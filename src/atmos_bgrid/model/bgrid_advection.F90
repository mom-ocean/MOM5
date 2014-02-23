
                   module bgrid_advection_mod

!-----------------------------------------------------------------------
!
!            performs vertical and horizontal advection
!              and negative tracer borrowing/filling
!
!-----------------------------------------------------------------------
!--------------------------- modules -----------------------------------
!-----------------------------------------------------------------------
use bgrid_horiz_mod       , only: horiz_grid_type
use bgrid_vert_mod        , only: vert_grid_type, compute_pres_depth
use bgrid_masks_mod       , only: grid_mask_type
use bgrid_prog_var_mod    , only: prog_var_type
use bgrid_polar_filter_mod, only: pfilt_control_type, polar_filter, &
                                  polar_filter_wind, TGRID
use bgrid_halo_mod        , only: update_halo, vel_flux_boundary,   &
                                  EAST, WEST, SOUTH, NORTH, NOPOLE, &
                                  TEMP, UWND, VWND, WIND, POLEONLY
use bgrid_change_grid_mod , only: change_grid, EQUAL, AREA, &
                                  UFLX_GRID, VFLX_GRID, TEMP_GRID, WIND_GRID
use vert_advection_mod    , only: vert_advection, vert_advection_end, &
                                  SECOND_CENTERED, FOURTH_CENTERED, &
                                  SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                  FINITE_VOLUME_LINEAR, FINITE_VOLUME_PARABOLIC, &
                                  FLUX_FORM, WEIGHTED_TENDENCY
!!!use horiz_advection_mod   , only: horiz_advection, &
!!!                                  FINITE_VOLUME_LINEAR_HORIZ=>FINITE_VOLUME_LINEAR, &
!!!                                  FINITE_VOLUME_PARABOLIC_HORIZ=>FINITE_VOLUME_PARABOLIC

use fms_mod, only: error_mesg, FATAL, write_version_number,      &
                   file_exist, open_namelist_file,               &
                   check_nml_error, close_file,                  &
                   mpp_pe, mpp_root_pe, stdlog, uppercase,       &
                   mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                   MPP_CLOCK_SYNC, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP
use mpp_mod, only: input_nml_file, mpp_max
use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_names, get_number_tracers
!-----------------------------------------------------------------------

implicit none
private

 public :: advection_init, advection, advection_end

!-----------------------------------------------------------------------
!------------ namelist: bgrid_advection_nml -------------

!  vert_advec_scheme_wind      The vertical advection scheme.
!  vert_advec_scheme_temp      Possible values are NONE, SECOND_CENTERED, FOURTH_CENTERED,
!  vert_advec_scheme_tracer    FINITE_VOLUME_LINEAR, or FINITE_VOLUME_PARABOLIC.
!                              Using finite volume schemes for momentum has not been
!                              tested and may produce poor results.

   character(len=24) ::  vert_advec_scheme_wind   = 'SECOND_CENTERED'
   character(len=24) ::  vert_advec_scheme_temp   = 'SECOND_CENTERED'
   character(len=24) ::  vert_advec_scheme_tracer = 'SECOND_CENTERED'

!  horiz_advec_scheme_wind     The horizontal advection scheme.
!  horiz_advec_scheme_temp     Possible values are NONE, SECOND_CENTERED, FOURTH_CENTERED.
!  horiz_advec_scheme_tracer

   character(len=24) :: horiz_advec_scheme_wind   = 'SECOND_CENTERED'
   character(len=24) :: horiz_advec_scheme_temp   = 'SECOND_CENTERED'
   character(len=24) :: horiz_advec_scheme_tracer = 'SECOND_CENTERED'

!  advec_weight_wind     Weights used for modified Euler-backward time differencing
!  advec_weight_temp     (i.e., when the scheme is SECOND_CENTERED or FOURTH_CENTERED).
!  advec_weight_tracer        0.0 = full Euler-forward (not recommended)
!                             1.0 = full Euler-backward

   real :: advec_weight_wind   = 0.7
   real :: advec_weight_temp   = 0.7
   real :: advec_weight_tracer = 0.7

!  num_fill_pass      The number of successive repetitions of the tracer borrowing scheme.
!                     This value applies to both the horizontal and vertical schemes.

   integer :: num_fill_pass = 1

!  temporary undocumented developer flags

logical :: compute_vert_wind_flux = .false.
character(len=16) ::  vert_vel_flux =  'area_weight'
character(len=16) :: horiz_vel_flux = 'equal_weight'

namelist /bgrid_advection_nml/ horiz_advec_scheme_wind,   vert_advec_scheme_wind,   &
                               horiz_advec_scheme_temp,   vert_advec_scheme_temp,   &
                               horiz_advec_scheme_tracer, vert_advec_scheme_tracer, &
                               advec_weight_wind, advec_weight_temp, advec_weight_tracer, &
                               num_fill_pass    ,&
                               compute_vert_wind_flux, vert_vel_flux, horiz_vel_flux

!-----------------------------------------------------------------------
!----- private data -----

! derived-type containing tracer filling parameters
 type trfill_control_type
    integer, pointer :: fill_scheme(:) =>NULL(), &
                        npass_horiz(:) =>NULL(), &
                        npass_vert(:)  =>NULL()
 end type trfill_control_type

! derived-type containing advection parameters
 type advec_control_type
    type(trfill_control_type) :: Fill
    integer, pointer  :: scheme(:,:) =>NULL()
    real   , pointer  :: weight(:) =>NULL()
    logical  ::  do_mask4_tmp, do_mask4_vel, &
                 do_finite_volume_tmp, do_finite_volume_trs
 end type advec_control_type

 integer, parameter :: NONE = 7000, &  ! may apply to advection or filling
            ! filling schemes
              EQUAL_BORROW  = 8005, & 
              BOTTOM_BORROW = 8006, &  ! not implemented
              GLOBAL_BORROW = 8007     ! not implemented

! advection control parameters
! set via values in namelist or field table
 type(advec_control_type),save :: Control

 real, parameter :: c4  = -1./6.       ! weights for 4th order centered schemes
 real, parameter :: c2  = 1. - 2.*c4

 logical :: do_log = .true.
 logical :: stability_check = .false.  ! perform stability check for horizontal
                                       ! finite volume (Lin-Rood) schemes

!  performance timing info
 integer :: id_tmp, id_trs, id_vel, id_fill
 logical :: do_clock_init = .true. 

 character(len=128) :: version='$Id: bgrid_advection.F90,v 19.0 2012/01/06 19:52:53 fms Exp $'
 character(len=128) :: tagname='$Name: tikal $'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine advection ( Pfilt, Hgrid, Vgrid, Masks,   &
                       pfilt_opt, dt, dpde_old, dpde,         &
                       few, fns, etadot, u, v, t, Var, Var_dt )

type(pfilt_control_type), intent(in) :: Pfilt
type(horiz_grid_type), intent(inout) :: Hgrid
type (vert_grid_type), intent(in)    :: Vgrid
type (grid_mask_type), intent(in)    :: Masks
              integer, intent(in)    :: pfilt_opt
                 real, intent(in)    :: dt
real, intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) ::  &
                                            dpde_old, dpde, u, v, t
real, intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: &
                                            few, fns, etadot
type (prog_var_type), intent(in)    :: Var
type (prog_var_type), intent(inout) :: Var_dt

!-----------------------------------------------------------------------
!
!  Pfilt    = polar filter constants
!  Hgrid    = horizontal grid constants
!  Vgrid    = vertical grid constants
!  Masks    = grid masking constants
!  pf_opt   = polar filter option flag
!  dt       = adjustment time step
!  dpde_old = pressure thickness of model layers at the end of the
!              last advective time step
!  dpde     = current pressure thickness of model layers
!  few, fns = zonal, meridional mass fluxes (Pa-m2/s)
!              (a summation since the last call to advection)
!  etadot   = vertical mass flux (Pa/s) (summation since the last call to advection)
!  u, v, t  = prognostic variables at the end of the last advective time
!              step, note that tracers have not been updated since the last
!              advective time step, therefore, r = Var%r + dt*Var_dt%r
!  Var      = prognostic variables at the end of the last dynamics time
!              step, current values would be, uc = Var%u + dt*Var_dt%u
!  Var_dt   = prognostic variable tendencies, accumulated since the
!              variable were updated in Var
!
!-----------------------------------------------------------------------
  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  Vgrid%nlev) :: dpdt, dpde_xy, r, uc, vc
  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  Vgrid%nlev,2) :: mask4
  integer :: i, j, k, n, is, ie, js, je
!-----------------------------------------------------------------------

!  ---- update halos for fluxes ----

   call update_halo (Hgrid, TEMP, few)
   call update_halo (Hgrid, VWND, fns)

!  ---- compute pressure tendency in each layer ----

   dpdt = (dpde - dpde_old) / dt

! ---- stability diagnostic (CFL for horizontal finite volume schemes) ----
   if (stability_check) call stability_diag ( Hgrid, dt, few, dpde_old, dpde )

!-----------------------------------------------------------------------
!------------------ advection of mass fields ---------------------------
!-----------------------------------------------------------------------

!------ initialize fourth-order mask ----

   if (Control%do_mask4_tmp) &
         call mask4_init (Hgrid, 1, Masks%sigma, Masks%Tmp%mask, mask4)

 ! compute wind components for horizontal finite-volume advection schemes
   if (Control%do_finite_volume_tmp .or. Control%do_finite_volume_trs) then
        call compute_advecting_wind (Hgrid, dpde_old, dpde, few, fns, uc, vc)
   else
        uc = 0.
        vc = 0.
   endif

!------ temperature advection ------

   call mpp_clock_begin (id_tmp)
   if (Masks%sigma) then
       call advect_mass (Hgrid, Pfilt, Control % scheme(:,0), &
                         Control % weight(0), pfilt_opt,      &
                         dt, dpdt, dpde, few, fns, etadot,    &
                         uc, vc, Var%t, t, Var_dt%t)
   else
       call advect_mass (Hgrid, Pfilt, Control % scheme(:,0), &
                         Control % weight(0), pfilt_opt,      &
                         dt, dpdt, dpde, few, fns, etadot,    &
                         uc, vc, Var%t, t, Var_dt%t,          &
                         mask=Masks%Tmp%mask, mask4=mask4)
   endif
   call mpp_clock_end   (id_tmp)

!------ tracer advection -----------
!  (need to pass all tracers to update pressure tendency part)

   call mpp_clock_begin (id_trs)
   do n = 1, Var_dt%ntrace

      r = Var%r(:,:,:,n) + dt*Var_dt%r(:,:,:,n)

      if (Masks%sigma) then
          call advect_mass (Hgrid, Pfilt, Control%scheme(:,n),   &
                            Control % weight(n),                 &
                            pfilt_opt, dt, dpdt, dpde,           &
                            few, fns, etadot, uc, vc,            &
                            Var%r(:,:,:,n), r, Var_dt%r(:,:,:,n) )
      else
          call advect_mass (Hgrid, Pfilt, Control%scheme(:,n),   &
                            Control % weight(n),                 &
                            pfilt_opt, dt, dpdt, dpde,           &
                            few, fns, etadot, uc, vc,            &
                            Var%r(:,:,:,n), r, Var_dt%r(:,:,:,n),&
                            mask=Masks%Tmp%mask, mask4=mask4     )
      endif

   enddo
   call mpp_clock_end   (id_trs)

!------ tracer hole filling -------
!  (remove negative tracer values with vert/horiz borrowing)

     call mpp_clock_begin (id_fill)
     call vert_borrow ( dt, dpde, Var%r(:,:,:,1:Var_dt%ntrace), &
                               Var_dt%r(:,:,:,1:Var_dt%ntrace), &
                                  iters = Control % Fill % npass_vert  )

     if (Masks%sigma) then
         call horiz_borrow ( Hgrid, dt, dpde,                 &
                                Var%r(:,:,:,1:Var_dt%ntrace), &
                             Var_dt%r(:,:,:,1:Var_dt%ntrace), &
                               iters = Control % Fill % npass_horiz  )
     else
         call horiz_borrow ( Hgrid, dt, dpde,                 &
                                Var%r(:,:,:,1:Var_dt%ntrace), &
                             Var_dt%r(:,:,:,1:Var_dt%ntrace), &
                                mask = Masks%Tmp%mask,        &
                               iters = Control % Fill % npass_horiz  )
     endif
     call mpp_clock_end   (id_fill)

!-----------------------------------------------------------------------
!------------------ advection of momentum fields -----------------------
!-----------------------------------------------------------------------

   call mpp_clock_begin (id_vel)

   is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
   js = Hgrid % Vel % js;  je = Hgrid % Vel % je

  !--- mass fluxes between velocity points ---

   if (trim(horiz_vel_flux(1:5)) == 'equal') then
      call change_grid (Hgrid, UFLX_GRID, VFLX_GRID, few, few, weight=EQUAL)
      call change_grid (Hgrid, VFLX_GRID, UFLX_GRID, fns, fns, weight=EQUAL)
   else
      call change_grid (Hgrid, UFLX_GRID, VFLX_GRID, few, few, weight=AREA)
      call change_grid (Hgrid, VFLX_GRID, UFLX_GRID, fns, fns, weight=AREA)
   endif

 ! compute vertical flux for momentum from horizonal momentum fluxes
   if (compute_vert_wind_flux) then
       call compute_etadot_vel (Hgrid, Vgrid, Masks, few, fns, etadot )
   endif

 ! no flux across pole
   call vel_flux_boundary (Hgrid, fns)

  !--- interpolate mass fields to momentum grid ---
   if (.not.compute_vert_wind_flux) then
       if (trim(vert_vel_flux(1:4)) == 'area') then
          call change_grid (Hgrid, TEMP_GRID, WIND_GRID, etadot, etadot, weight=AREA)
       else
          call change_grid (Hgrid, TEMP_GRID, WIND_GRID, etadot, etadot, weight=EQUAL)
       endif
   endif
   ! always area-weighted (the default)
   call change_grid (Hgrid, TEMP_GRID, WIND_GRID, dpde, dpde_xy)
   call change_grid (Hgrid, TEMP_GRID, WIND_GRID, dpdt,    dpdt)


 ! advection of momentum
 ! determine whether step-mountain mask are needed
   if (Control%do_mask4_vel .and. .not.Masks%sigma) then
       call mask4_init (Hgrid, 2, Masks%sigma, Masks%Vel%mask, mask4)

       call advect_vel (Hgrid, Pfilt, Control % scheme(:,-1), &
                        Control % weight(-1),                 &
                        pfilt_opt, dt, dpdt, dpde_xy,         &
                        few, fns, etadot, Var%u, Var%v, u, v, &
                        Var_dt%u, Var_dt%v,                   &
                        mask=Masks%Vel%mask, mask4=mask4)
   else
       ! sigma case
       call advect_vel (Hgrid, Pfilt, Control % scheme(:,-1), &
                        Control % weight(-1),                 &
                        pfilt_opt, dt, dpdt, dpde_xy,         &
                        few, fns, etadot, Var%u, Var%v, u, v, &
                        Var_dt%u, Var_dt%v)
   endif

   call mpp_clock_end (id_vel)
!-----------------------------------------------------------------------
!------ done with fluxes - zero-out ??? -------

   few    = 0.0
   fns    = 0.0
   etadot = 0.0

!-----------------------------------------------------------------------

end subroutine advection

!#######################################################################

subroutine advect_mass ( Hgrid, Pfilt, scheme, coeff, fopt, dt, &
                         dpdt, dpn, few, fns, fud, uc, vc,      &
                         ro, r, r_dt, mask, mask4 )

  type(horiz_grid_type),    intent(inout) :: Hgrid
  type(pfilt_control_type), intent(in) :: Pfilt
  integer,                  intent(in) :: scheme(2)
  integer,                  intent(in) :: fopt
  real, intent(in)                     :: coeff, dt
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpdt, dpn,&
                                           few, fns, fud, ro, r, uc, vc
  real, intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: r_dt
  real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:),   optional :: mask
  real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:), optional :: mask4

!-----------------------------------------------------------------------

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  size(r,3)) :: rst, rst_dt, rst_dt_v, rdpdt

  real    :: rcoef
  integer :: pass, npass, j, k
  integer :: vert_scheme, horiz_scheme
  integer, parameter :: FINITE_VOLUME = 1234567 !can not match other numbers?

!-----------------------------------------------------------------------

 ! set the vertical and horizontal differencing scheme
    horiz_scheme = scheme(1)
     vert_scheme = scheme(2)
 ! rename horizontal finite volume schemes to correct horizontal name
 !  if (horiz_scheme == FINITE_VOLUME_LINEAR)    horiz_scheme = FINITE_VOLUME_LINEAR_HORIZ
 !  if (horiz_scheme == FINITE_VOLUME_PARABOLIC) horiz_scheme = FINITE_VOLUME_PARABOLIC_HORIZ

!-----------------------------------------------------------------------

   rst = ro + dt * r_dt  ! use updated value
   rdpdt = r*dpdt
   rcoef = coeff

                      npass = 1
   if (rcoef > 1.e-4) npass = 2
 ! no advection, psdt correction only
   if (vert_scheme ==  NONE) then
       rst_dt_v = 0.0
       if (horiz_scheme == NONE) npass = 1
   endif
 ! both vert and horiz are using finite volume
   if ( (vert_scheme == FINITE_VOLUME_LINEAR .or.       &
         vert_scheme == FINITE_VOLUME_PARABOLIC) .and.  &
       (horiz_scheme == FINITE_VOLUME_LINEAR .or. &
        horiz_scheme == FINITE_VOLUME_PARABOLIC) ) then
      !(horiz_scheme == FINITE_VOLUME_LINEAR_HORIZ .or. &
      ! horiz_scheme == FINITE_VOLUME_PARABOLIC_HORIZ) ) then
          npass = 2
          rcoef = 1.
   endif

   do pass = 1, npass

     !---- vertical tendency ----
         ! finite volume schemes on pass 1 only
      if (vert_scheme == SECOND_CENTERED .or. vert_scheme == SECOND_CENTERED_WTS .or. &
          vert_scheme == FOURTH_CENTERED .or. vert_scheme == FOURTH_CENTERED_WTS .or. &
        ((vert_scheme == FINITE_VOLUME_LINEAR .or. &
          vert_scheme == FINITE_VOLUME_PARABOLIC) .and. pass == 1)) then
         call vert_advection ( dt, fud, dpn, rst, rst_dt_v,   &
                               mask=mask, scheme=vert_scheme, &
                               form=FLUX_FORM, flags=WEIGHTED_TENDENCY )
      endif

     !---- horizontal tendency ----

      if (horiz_scheme == SECOND_CENTERED) then
         call advect_mass_horiz_2nd (Hgrid, few, fns, rst, rst_dt)
         if ( fopt >= 1 ) call polar_filter (Pfilt, rst_dt, TGRID, mask)
         call update_halo (Hgrid, TEMP, rst_dt)

      else if (horiz_scheme == FOURTH_CENTERED) then
         call advect_mass_horiz_4th (Hgrid, few, fns, rst, rst_dt, mask4)
         if ( fopt >= 1 ) call polar_filter (Pfilt, rst_dt, TGRID, mask)
         call update_halo (Hgrid, TEMP, rst_dt)

    ! compute horizontal finite volume scheme on 2nd pass
      else if ((horiz_scheme == FINITE_VOLUME_LINEAR .or. &
                horiz_scheme == FINITE_VOLUME_PARABOLIC) .and. pass == 2) then
        ! temporary error check
         call error_mesg ('bgrid_advection_mod', &
                          'horizontal finite volume schemes &
                          &are not implemented with this release', FATAL)
        !call horiz_advection ( Hgrid%Tmp%Domain, dt, Hgrid%Tmp%dx, Hgrid%Tmp%dy, &
        !                       uc, vc, few, fns, rst, rst_dt, scheme=horiz_scheme )
       ! halos updated by finite volume scheme
       ! need update only poles (model b.c. differs)
         call update_halo (Hgrid, TEMP, rst_dt, flags=POLEONLY)

      else
            rst_dt = 0.0
      endif
      

     !---- combine vertical and horizontal tendencies ----
     !---- adjust for pressure tendency ----
        rst_dt = (rst_dt + rst_dt_v - rdpdt) / dpn
      ! apply step-mountain mask?
        if (present(mask)) rst_dt = rst_dt*mask

     !---- compute new value or return tendency ----
        if (pass < npass) then
            rst = ro + dt * (r_dt + rcoef*rst_dt)
        else
            r_dt = r_dt + rst_dt
        endif

   enddo

!-----------------------------------------------------------------------

end subroutine advect_mass

!#######################################################################

subroutine advect_vel ( Hgrid, Pfilt, scheme, coeff, fopt, dt,  &
                        dpdt, dpn, few, fns, fud,               &
                        uo, vo, u, v, u_dt, v_dt, mask, mask4  )

  type(horiz_grid_type),    intent(inout) :: Hgrid
  type(pfilt_control_type), intent(in) :: Pfilt
  integer,                  intent(in) :: scheme(2)
  integer,                  intent(in) :: fopt
  real, intent(in)                     :: coeff, dt
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpdt, dpn,&
                                           few, fns, fud, uo, vo, u, v
  real, intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u_dt, v_dt
  real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:),   optional :: mask
  real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:), optional :: mask4

!-----------------------------------------------------------------------

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  size(u,3)) :: ust_dt, ust_dt_v, &
                                vst_dt, vst_dt_v
! store ust & vst components together to create larger halo updates
  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  size(u,3),2) :: uvst
  integer, parameter :: UCOMP=1, VCOMP=2

  integer :: vert_scheme, horiz_scheme
  integer :: pass, npass, i, j, k, is, ie, js, je, nlev

!-----------------------------------------------------------------------

 ! set the vertical differencing scheme
   horiz_scheme = scheme(1)
    vert_scheme = scheme(2)

!-----------------------------------------------------------------------

   is = Hgrid%Vel%is; ie = Hgrid%Vel%ie
   js = Hgrid%Vel%js; je = Hgrid%Vel%je
   nlev = size(u,3)

   uvst(:,:,:,UCOMP) = uo + dt * u_dt   ! use updated values
   uvst(:,:,:,VCOMP) = vo + dt * v_dt

 ! tendency halos are not up-to-date
 ! need to update halos
   call update_halo ( Hgrid, WIND, uvst )

                      npass = 1
   if (coeff > 1.e-4) npass = 2

   do pass = 1, npass

     !---- vertical tendency ----
         ! finite volume schemes on pass 1 only
      if (vert_scheme == SECOND_CENTERED .or. vert_scheme == SECOND_CENTERED_WTS .or. &
          vert_scheme == FOURTH_CENTERED .or. vert_scheme == FOURTH_CENTERED_WTS .or. &
        ((vert_scheme == FINITE_VOLUME_LINEAR .or. &
          vert_scheme == FINITE_VOLUME_PARABOLIC) .and. pass == 1)) then
          call vert_advection ( dt, fud, dpn, uvst(:,:,:,UCOMP), ust_dt_v, &
                                mask=mask, scheme=vert_scheme,             &
                                form=FLUX_FORM, flags=WEIGHTED_TENDENCY )
          call vert_advection ( dt, fud, dpn, uvst(:,:,:,VCOMP), vst_dt_v, &
                                mask=mask, scheme=vert_scheme,             &
                                form=FLUX_FORM, flags=WEIGHTED_TENDENCY )
      else if (vert_scheme == NONE .and. pass == 1) then
          ust_dt_v = 0.0
          vst_dt_v = 0.0
      endif

     !---- horizontal tendency ----
      if (horiz_scheme == SECOND_CENTERED) then
        call advect_vel_horiz_2nd (Hgrid, few, fns, uvst(:,:,:,UCOMP), &
                                   uvst(:,:,:,VCOMP), ust_dt, vst_dt   )
      
      else if (horiz_scheme == FOURTH_CENTERED) then
        call advect_vel_horiz_4th (Hgrid, few, fns,         &
                               uvst(:,:,:,UCOMP), uvst(:,:,:,VCOMP), &
                               ust_dt, vst_dt, mask4 )
      endif
      
     ! polar filter horizontal tendency
        if (fopt >= 2) then
           call polar_filter_wind (Pfilt, ust_dt, vst_dt, mask)
        endif

     !---- combine vertical and horizontal tendencies ----
     !---- adjust for pressure tendency ----
      do k = 1, nlev
      do j = js, je
      do i = is, ie
         ust_dt(i,j,k) = (ust_dt(i,j,k) + ust_dt_v(i,j,k) - &
                       u(i,j,k)*dpdt(i,j,k)) / dpn(i,j,k)
         vst_dt(i,j,k) = (vst_dt(i,j,k) + vst_dt_v(i,j,k) - &
                       v(i,j,k)*dpdt(i,j,k)) / dpn(i,j,k)
      enddo
      enddo
      enddo
    ! apply step-mountain mask?
      if (present(mask)) then
        do k = 1, nlev
        ust_dt(is:ie,js:je,k) = ust_dt(is:ie,js:je,k)*mask(is:ie,js:je,k)
        vst_dt(is:ie,js:je,k) = vst_dt(is:ie,js:je,k)*mask(is:ie,js:je,k)
        enddo
      endif

     !---- compute new value or return tendency ----
      if (pass < npass) then
         do k = 1, nlev
         do j = js, je
         do i = is, ie
           uvst(i,j,k,UCOMP) = uo(i,j,k) + dt * (u_dt(i,j,k) + coeff*ust_dt(i,j,k))
           uvst(i,j,k,VCOMP) = vo(i,j,k) + dt * (v_dt(i,j,k) + coeff*vst_dt(i,j,k))
         enddo
         enddo
         enddo
         call update_halo (Hgrid, WIND, uvst)
      else
         do k = 1, nlev
         do j = js, je
         do i = is, ie
           u_dt(i,j,k) = u_dt(i,j,k) + ust_dt(i,j,k)
           v_dt(i,j,k) = v_dt(i,j,k) + vst_dt(i,j,k)
         enddo
         enddo
         enddo
       ! NOTE: halos not updated for output tendencies
       !  do this in the dynamical core for efficiency
      endif
   
   enddo

!-----------------------------------------------------------------------

end subroutine advect_vel

!#######################################################################

 subroutine advect_mass_horiz_2nd (Hgrid, few, fns, r, rdt)

  type(horiz_grid_type),    intent(in) :: Hgrid
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, r
  real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: rdt

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: frew, frns
  integer :: i, j, k

! second order centered differencing on B-grid temperature grid
! constant grid spacing assumed
! minimum halo size = 1

      do k = 1, size(r,3)

       !--- horizontal fluxes ---
        do j = Hgrid%Tmp%js-1, Hgrid%Tmp%je
        do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
          frew(i,j) = few(i,j,k) * (r(i+1,j,k) + r(i,j,k))
          frns(i,j) = fns(i,j,k) * (r(i,j+1,k) + r(i,j,k))
        enddo
        enddo

       !--- horizontal advective tendency ---
        do j = Hgrid%Tmp%js, Hgrid%Tmp%je
        do i = Hgrid%Tmp%is, Hgrid%Tmp%ie
           rdt(i,j,k) = -0.5*Hgrid%Tmp%rarea(j)*       &
                           ((frew(i  ,j)+frns(i,j  ))  &
                           -(frew(i-1,j)+frns(i,j-1)))
        enddo
        enddo

      enddo

 end subroutine advect_mass_horiz_2nd

!#######################################################################

 subroutine advect_mass_horiz_4th (Hgrid, few, fns, r, rdt, mask4)

  type(horiz_grid_type),    intent(inout) :: Hgrid
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, r
  real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: rdt
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:,:), optional :: &
                                                           mask4

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, &
                  size(r,3)) :: rew, rns, frew, frns
  real, dimension(Hgrid%ilb:Hgrid%iub) :: x2, x4
  real, dimension(Hgrid%ilb:Hgrid%iub) :: y2, y4
  integer :: i, j, k, is, ie, js, je

! fourth order centered differencing on B-grid temperature grid
! constant grid spacing assumed
! assumed halo size = 1

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

      rns(:,Hgrid%jub,:) = 0.0
      do k = 1, size(r,3)
        do j = Hgrid%jlb,Hgrid%jub
        do i = Hgrid%ilb,Hgrid%iub-1
           rew(i,j,k) = r(i+1,j,k)+r(i,j,k)
        enddo
        enddo
        do j = Hgrid%jlb,Hgrid%jub-1
        do i = Hgrid%ilb,Hgrid%iub
           rns(i,j,k) = r(i,j+1,k)+r(i,j,k)
        enddo
        enddo
      enddo
      ! assumed halosize = 1
      call update_halo (Hgrid,TEMP,rew,halos=EAST)
      call update_halo (Hgrid,UWND,rns,halos=NORTH+NOPOLE)

!     ---- horizontal fluxes ----

      if (present(mask4)) then
         frew=0.0;  frns=0.0
         do k = 1, size(r,3)
           do j = js,          je
           do i = Hgrid%ilb+1, ie
             x4(i) = mask4(i,j,k,1)
             x2(i) = 1.0 - 2.*x4(i)
             frew(i,j,k) = few(i,j,k) * ( x2(i) * rew(i,j,k) &
                          + x4(i) * (rew(i-1,j,k) + rew(i+1,j,k)) )
           enddo
           enddo
           do j = Hgrid%jlb+1, je
           do i = is,          ie
             y4(i) = mask4(i,j,k,2)
             y2(i) = 1.0 - 2.*y4(i)
             frns(i,j,k) = fns(i,j,k) * ( y2(i) * rns(i,j,k) &
                          + y4(i) * (rns(i,j-1,k) + rns(i,j+1,k)) )
           enddo
           enddo
         enddo
      else
         frew=0.0;  frns=0.0
         do k = 1, size(r,3)
           do j = js,          je
           do i = Hgrid%ilb+1, ie
             frew(i,j,k) = few(i,j,k) * ( c2 * rew(i,j,k) &
                          + c4 * (rew(i-1,j,k) + rew(i+1,j,k)) )
           enddo
           enddo
           do j = Hgrid%jlb+1, je
           do i = is,          ie
             frns(i,j,k) = fns(i,j,k) * ( c2 * rns(i,j,k) &
                          + c4 * (rns(i,j-1,k) + rns(i,j+1,k)) )
           enddo
           enddo
         enddo
      endif

      ! assumed halosize = 1
      call update_halo (Hgrid,TEMP,frew,halos=WEST)
      call update_halo (Hgrid,VWND,frns,halos=SOUTH)

!     ---- horizontal advective tendency ----

      do k = 1, size(r,3)
      do j = js, je
      do i = is, ie
         rdt(i,j,k) = -0.5*Hgrid%Tmp%rarea(j)*           &
                         ((frew(i  ,j,k)+frns(i,j  ,k))  &
                         -(frew(i-1,j,k)+frns(i,j-1,k)))
      enddo
      enddo
      enddo


 end subroutine advect_mass_horiz_4th

!#######################################################################

 subroutine advect_vel_horiz_2nd (Hgrid, few, fns, u, v, udt, vdt)

  type(horiz_grid_type),    intent(in) :: Hgrid
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, &
                                                            u, v
  real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: udt, vdt

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
                                             fuew, funs, fvew, fvns
  integer :: i, j, k

! second order centered differencing on B-grid momentum grid
! constant grid spacing assumed
! minimum halo size = 1

     do k = 1, size(u,3)

          !--- horizontal fluxes ---
           do j = Hgrid%Vel%js, Hgrid%Vel%je
           do i = Hgrid%Vel%is, Hgrid%Vel%ie+1
             fuew(i,j) = few(i,j,k) * (u(i-1,j,k) + u(i,j,k))
             fvew(i,j) = few(i,j,k) * (v(i-1,j,k) + v(i,j,k))
           enddo
           enddo
           do j = Hgrid%Vel%js, Hgrid%Vel%je+1
           do i = Hgrid%Vel%is, Hgrid%Vel%ie
             funs(i,j) = fns(i,j,k) * (u(i,j-1,k) + u(i,j,k))
             fvns(i,j) = fns(i,j,k) * (v(i,j-1,k) + v(i,j,k))
           enddo
           enddo

          !--- horizontal advective tendency ----
           do j = Hgrid%Vel%js, Hgrid%Vel%je
           do i = Hgrid%Vel%is, Hgrid%Vel%ie
              udt(i,j,k) = -0.5*Hgrid%Vel%rarea(j)*       &
                              ((fuew(i+1,j)+funs(i,j+1))  &
                              -(fuew(i  ,j)+funs(i,j  )))
              vdt(i,j,k) = -0.5*Hgrid%Vel%rarea(j)*       &
                              ((fvew(i+1,j)+fvns(i,j+1))  &
                              -(fvew(i  ,j)+fvns(i,j  )))
           enddo
           enddo

     enddo

 end subroutine advect_vel_horiz_2nd

!#######################################################################

 subroutine advect_vel_horiz_4th (Hgrid, few, fns, u, v, udt, vdt, &
                                  mask4)

  type(horiz_grid_type),    intent(inout) :: Hgrid
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, u, v
  real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: udt, vdt
  real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:,:), optional :: &
                                                           mask4

! compress u & v components into 4th dimension
! this will create larger and more efficient halo updates
  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, &
                  size(u,3),2) :: uvew, uvns, fuvew, fuvns
  integer, parameter :: UCOMP=1, VCOMP=2

  real, dimension(Hgrid%ilb:Hgrid%iub) :: x2, x4
  real, dimension(Hgrid%ilb:Hgrid%iub) :: y2, y4
  real :: z2, z4
  integer :: i, j, k, is, ie, js, je

! fourth order centered differencing on B-grid momentum grid
! constant grid spacing assumed
! assumed halo size = 1

   is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
   js = Hgrid % Vel % js;  je = Hgrid % Vel % je


      do k = 1, size(u,3)
        do j = Hgrid%jlb  ,Hgrid%jub
        do i = Hgrid%ilb+1,Hgrid%iub
           uvew(i,j,k,UCOMP) = u(i-1,j,k)+u(i,j,k)
           uvew(i,j,k,VCOMP) = v(i-1,j,k)+v(i,j,k)
        enddo
        enddo
        do j = Hgrid%jlb+1,Hgrid%jub
        do i = Hgrid%ilb  ,Hgrid%iub
           uvns(i,j,k,UCOMP) = u(i,j-1,k)+u(i,j,k)
           uvns(i,j,k,VCOMP) = v(i,j-1,k)+v(i,j,k)
        enddo
        enddo
      enddo
      ! assumed halosize = 1
      call update_halo (Hgrid,WIND,uvew,halos=WEST)
      call update_halo (Hgrid,TEMP,uvns,halos=SOUTH)


     !---- horizontal fluxes ----

      if (present(mask4)) then
         do k = 1, size(u,3)
           do j = js, je
           do i = is, Hgrid%iub-1
             x4(i) = mask4(i,j,k,1)
             x2(i) = 1.0 - 2.*x4(i)
             fuvew(i,j,k,UCOMP) = few(i,j,k) * ( x2(i) * uvew(i,j,k,UCOMP) &
                            + x4(i) * (uvew(i-1,j,k,UCOMP) + uvew(i+1,j,k,UCOMP)) )
             fuvew(i,j,k,VCOMP) = few(i,j,k) * ( x2(i) * uvew(i,j,k,VCOMP) &
                            + x4(i) * (uvew(i-1,j,k,VCOMP) + uvew(i+1,j,k,VCOMP)) )
           enddo
           enddo
           do j = js, Hgrid%jub-1
           do i = is, ie
             y4(i) = mask4(i,j,k,2)
             y2(i) = 1.0 - 2.*y4(i)
             fuvns(i,j,k,UCOMP) = fns(i,j,k) * ( y2(i) * uvns(i,j,k,UCOMP) &
                                + y4(i) * (uvns(i,j-1,k,UCOMP) + uvns(i,j+1,k,UCOMP)) )
             fuvns(i,j,k,VCOMP) = fns(i,j,k) * ( y2(i) * uvns(i,j,k,VCOMP) &
                                + y4(i) * (uvns(i,j-1,k,VCOMP) + uvns(i,j+1,k,VCOMP)) )
           enddo
           enddo
         enddo
      else
         do k = 1, size(u,3)
           do j = js, je
           do i = is, Hgrid%iub-1
             fuvew(i,j,k,UCOMP) = few(i,j,k) * ( c2 * uvew(i,j,k,UCOMP) &
                            + c4 * (uvew(i-1,j,k,UCOMP) + uvew(i+1,j,k,UCOMP)) )
             fuvew(i,j,k,VCOMP) = few(i,j,k) * ( c2 * uvew(i,j,k,VCOMP) &
                            + c4 * (uvew(i-1,j,k,VCOMP) + uvew(i+1,j,k,VCOMP)) )
           enddo
           enddo
           do j = js, Hgrid%jub-1
           ! second order near poles
             z4 = c4
             if (j <= Hgrid%Vel%jsg+1) z4 = 0.
             if (j >= Hgrid%Vel%jeg  ) z4 = 0.
             z2 = 1.-2.*z4
           do i = is, ie
             fuvns(i,j,k,UCOMP) = fns(i,j,k) * ( z2 * uvns(i,j,k,UCOMP) &
                                + z4 * (uvns(i,j-1,k,UCOMP) + uvns(i,j+1,k,UCOMP)) )
             fuvns(i,j,k,VCOMP) = fns(i,j,k) * ( z2 * uvns(i,j,k,VCOMP) &
                                + z4 * (uvns(i,j-1,k,VCOMP) + uvns(i,j+1,k,VCOMP)) )
           enddo
           enddo
         enddo
      endif
      ! assumed halosize = 1
      call update_halo (Hgrid,WIND,fuvew,halos=EAST)
      call update_halo (Hgrid,TEMP,fuvns,halos=NORTH)

     !---- horizontal advective tendency ----

      do k = 1, size(u,3)
      do j = js, je
      do i = is, ie
         udt(i,j,k) = -0.5*Hgrid%Vel%rarea(j)*                        &
                        ((fuvew(i+1,j,k,UCOMP)+fuvns(i,j+1,k,UCOMP))  &
                        -(fuvew(i  ,j,k,UCOMP)+fuvns(i,j  ,k,UCOMP)))
         vdt(i,j,k) = -0.5*Hgrid%Vel%rarea(j)*                        &
                        ((fuvew(i+1,j,k,VCOMP)+fuvns(i,j+1,k,VCOMP))  &
                        -(fuvew(i  ,j,k,VCOMP)+fuvns(i,j  ,k,VCOMP)))
      enddo
      enddo
      enddo


 end subroutine advect_vel_horiz_4th

!#######################################################################
! computes advecting wind for horizontal finite-volume advection

 subroutine compute_advecting_wind ( Hgrid, dpo, dpn, few, fns, uc, vc )
 type(horiz_grid_type), intent(inout) :: Hgrid
 real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dpo, dpn, few, fns
 real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: uc, vc

 real :: dp
 integer :: i, j, k

   do k = 1, size(uc,3)
   do j = Hgrid%jlb, Hgrid%jub
   do i = Hgrid%ilb, Hgrid%iub-1
     dp = (dpo(i,j,k)+dpn(i,j,k)+dpo(i+1,j,k)+dpn(i+1,j,k))*0.25
     uc(i,j,k) = few(i,j,k)/(dp*Hgrid%Tmp%dy)
   enddo
   enddo
   enddo

   do k = 1, size(uc,3)
   do j = Hgrid%jlb, Hgrid%jub-1
   do i = Hgrid%ilb, Hgrid%iub
     dp = (dpo(i,j,k)+dpn(i,j,k)+dpo(i,j+1,k)+dpn(i,j+1,k))*0.25
     vc(i,j,k) = fns(i,j,k)/(dp*Hgrid%Vel%dx(j))
   enddo
   enddo
   enddo

   call update_halo ( Hgrid, TEMP, uc, halos=EAST )
   call update_halo ( Hgrid, VWND, vc, halos=NORTH )

 end subroutine compute_advecting_wind

!#######################################################################
!--- stability diagnostic (CFL for horizontal finite volume scheme) ----

 subroutine stability_diag ( Hgrid, dt, few, dpo, dpn )
 type(horiz_grid_type), intent(inout) :: Hgrid
 real, intent(in)                     :: dt
 real, intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, dpo, dpn

 real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: uc
 real    :: cflmax, cfl, dpa
 integer :: i, j, k

    cflmax = 0.
    do k = 1, size(few,3)
     ! compute zonal wind
       do j = Hgrid%Tmp%js  , Hgrid%Tmp%je
       do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
         dpa = (dpo(i,j,k)+dpn(i,j,k)+dpo(i+1,j,k)+dpn(i+1,j,k))*0.25
         uc(i,j) = few(i,j,k)/(dpa*Hgrid%Tmp%dy)
       enddo     
       enddo   
       do j = Hgrid%Tmp%js, Hgrid%Tmp%je
       do i = Hgrid%Tmp%is, Hgrid%Tmp%ie
          cfl = abs(uc(i,j)-uc(i-1,j))*dt/Hgrid%Tmp%dx(j)
          cflmax = max(cflmax,cfl)
       enddo      
       enddo   
    enddo   
    call mpp_max (cflmax)
    if (mpp_pe() == mpp_root_pe()) then
        if (cflmax > 1.0) then
            print *, 'x-axis stability violated, cfl = ', cflmax
        endif       
    endif       

 end subroutine stability_diag

!#######################################################################

 subroutine mask4_init (Hgrid, grid, sigma, mask, mask4)

 type(horiz_grid_type), intent(inout) :: Hgrid
 integer, intent(in)                                      :: grid
 logical, intent(in)                                      :: sigma
 real, intent(in),   dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: mask
 real, intent(out),  dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: mask4
 real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                 size(mask,3)) :: mew, mns

 integer :: i,j,k,is,ie,js,je,n


! horizontal masks

  select case (grid)

  case (1)
     if (sigma) then
        mask4(:,:,:,1:2) = c4
     else

      ! initialize
        mask4 = 0.0
        mns(:,Hgrid%jub,:) = 0.0

        do k = 1, size(mask,3)
          do j = Hgrid%jlb,Hgrid%jub
          do i = Hgrid%ilb,Hgrid%iub-1
             mew(i,j,k) = mask(i+1,j,k)*mask(i,j,k)
          enddo
          enddo
          do j = Hgrid%jlb,Hgrid%jub-1
          do i = Hgrid%ilb,Hgrid%iub
             mns(i,j,k) = mask(i,j+1,k)*mask(i,j,k)
          enddo
          enddo
        enddo
        ! assumed halosize = 1
        call update_halo (Hgrid,TEMP,mew,halos=EAST)
        call update_halo (Hgrid,UWND,mns,halos=NORTH,flags=NOPOLE)

         do k = 1, size(mask,3)
         ! x-axis
           do j = Hgrid%Tmp%js, Hgrid%Tmp%je
           do i = Hgrid%ilb+1,  Hgrid%Tmp%ie
             !mask4(i,j,k,1) = mew(i-1,j,k)<.01 .or. mew(i,j,k)<.01 .or. mew(i+1,j,k)<.01
             !mask4(i,j,k,1) = mew(i-1,j,k) * mew(i,j,k) * mew(i+1,j,k)
              mask4(i,j,k,1) = c4 * mew(i-1,j,k) * mew(i+1,j,k)
           enddo
           enddo
         ! y-axis
           do j = Hgrid%jlb+1,  Hgrid%Tmp%je
           do i = Hgrid%Tmp%is, Hgrid%Tmp%ie
             !mask4(i,j,k,2) = mns(i,j-1,k)<.01 .or. mns(i,j,k)<.01 .or. mns(i,j-1,k)<.01
             !mask4(i,j,k,2) = mns(i,j-1,k) * mns(i,j,k) * mns(i,j+1,k)
              mask4(i,j,k,2) = c4 * mns(i,j-1,k) * mns(i,j+1,k)
           enddo
           enddo
         enddo
     endif

  case (2)
     if (sigma) then
        mask4(:,:,:,1:2) = c4
        mns = 1.0
     else
      ! initialize
        mask4 = 0.0
        do k = 1, size(mask,3)
          do j = Hgrid%jlb  ,Hgrid%jub
          do i = Hgrid%ilb+1,Hgrid%iub
             mew(i,j,k) = mask(i-1,j,k)*mask(i,j,k)
          enddo
          enddo
          do j = Hgrid%jlb+1,Hgrid%jub
          do i = Hgrid%ilb  ,Hgrid%iub
             mns(i,j,k) = mask(i,j-1,k)*mask(i,j,k)
          enddo
          enddo
        enddo
        ! assumed halosize = 1
        call update_halo (Hgrid,UWND,mew,halos=WEST)
        call update_halo (Hgrid,TEMP,mns,halos=SOUTH,flags=NOPOLE)
      endif
    ! use second order in meridional direction near poles
      call vel_flux_boundary (Hgrid,mns)

      if (.not.sigma) then
        do k = 1, size(mask,3)
          do j = Hgrid%Vel%js, Hgrid%Vel%je
          do i = Hgrid%Vel%is, Hgrid%iub-1
             mask4(i,j,k,1) = c4 * mew(i-1,j,k) * mew(i+1,j,k)
          enddo
          enddo
        enddo
      endif

      do k = 1, size(mask,3)
        do j = Hgrid%Vel%js, Hgrid%jub-1
        do i = Hgrid%Vel%is, Hgrid%Vel%ie
           mask4(i,j,k,2) = c4 * mns(i,j-1,k) * mns(i,j+1,k)
        enddo
        enddo
      enddo

  end select


 end subroutine mask4_init

!#######################################################################
! initialization routine

  subroutine advection_init ( Hgrid )

! INPUT:  Hgrid = horizontal grid constants
! RESULT: reads bgrid_advection_nml namelist
!         defines advection control parameters (stored in advec_control_type)

  type(horiz_grid_type), intent(in) :: Hgrid

 real    :: wt
 integer :: n, m, np, ntrace, unit, ierr, io, logunit
 character(len=32)  :: advec_methods(2) = (/ 'advec_horiz', 'advec_vert ' /)
 character(len=128) :: scheme, params, name
!-----------------------------------------------------------------------
! read namelist
   if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=bgrid_advection_nml, iostat=io)
       ierr = check_nml_error(io,'bgrid_advection_nml')
#else
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read (unit, nml=bgrid_advection_nml, iostat=io, end=5)
          ierr = check_nml_error (io, 'bgrid_advection_nml')
       enddo
 5     call close_file (unit)
#endif
   endif
   logunit = stdlog()
! write version, namelist info to log file
   if (do_log) then
      call write_version_number (version,tagname)
      if (mpp_pe() == mpp_root_pe()) write (logunit, nml=bgrid_advection_nml)
      do_log = .false.
   endif

 ! determine the number of prognostic tracers
   call get_number_tracers ( MODEL_ATMOS, num_prog=ntrace )

 ! allocate space for the control parameters
   allocate ( Control%scheme (2,-1:ntrace), &
              Control%weight   (-1:ntrace)  )
   allocate ( Control%Fill%fill_scheme (ntrace), &
              Control%Fill%npass_horiz (ntrace), &
              Control%Fill%npass_vert  (ntrace)  )

 ! set namelist values of control parameters

   Control%scheme(1,-1) = set_advec_scheme( horiz_advec_scheme_wind )
   Control%scheme(2,-1) = set_advec_scheme(  vert_advec_scheme_wind )
   Control%weight  (-1) = advec_weight_wind
   Control%scheme(1, 0) = set_advec_scheme( horiz_advec_scheme_temp )
   Control%scheme(2, 0) = set_advec_scheme(  vert_advec_scheme_temp )
   Control%weight  ( 0) = advec_weight_temp
 do n = 1, ntrace
   Control%scheme(1, n) = set_advec_scheme( horiz_advec_scheme_tracer )
   Control%scheme(2, n) = set_advec_scheme(  vert_advec_scheme_tracer )
   Control%weight  ( n) = advec_weight_tracer
   ! parameters for tracer filling/borrowing
   Control%Fill%npass_horiz(n) = num_fill_pass
   Control%Fill%npass_vert (n) = num_fill_pass
   Control%Fill%fill_scheme(n) = NONE
 enddo
   ! currently only one possible filling scheme
   if (maxval(Control%Fill%npass_horiz) > 0 .or. maxval(Control%Fill%npass_vert) > 0) &
                            Control%Fill%fill_scheme = EQUAL_BORROW

 ! process tracer table information for advection methods
   do n = 1, ntrace
      do m = 1, 2
         if (query_method(trim(advec_methods(m)), MODEL_ATMOS, n, scheme, params)) then
             Control%scheme(m,n) = set_advec_scheme( scheme )
             ! parse e-b weight
             if (Control%scheme(m,n) == SECOND_CENTERED .or. &
                 Control%scheme(m,n) == FOURTH_CENTERED) then
                   if (parse(params,'wt', wt) == 1) Control%weight(n) = wt
             endif
         endif
      enddo
   enddo

 ! error check for all Euler weights
   do n = -1, ntrace
        if (Control%weight(n) < 0.0 .or. Control%weight(n) > 1.0) &
                          call error_mesg ('bgrid_advection_mod', &
                            'E-B weight out of range [0,1]', FATAL)
   enddo

 ! process tracer table information for filling method
   do n = 1, ntrace
      if (query_method('filling', MODEL_ATMOS, n, scheme, params)) then
          if (uppercase(trim(scheme)) == 'LOCAL') then
              Control%Fill%fill_scheme(n) = EQUAL_BORROW
              Control%Fill%npass_horiz(n) = 1
              Control%Fill%npass_vert (n) = 1
         !else if (uppercase(trim(scheme)) == 'GLOBAL') then
         !    Control%Fill%fill_scheme(n) = GLOBAL_BORROW
          else if (uppercase(trim(scheme)) == 'NONE') then
              Control%Fill%fill_scheme(n) = NONE
              Control%Fill%npass_horiz(n) = 0
              Control%Fill%npass_vert (n) = 0
          else
              call error_mesg ('bgrid_advection_mod',  &
                   'invalid filling scheme, '//uppercase(trim(scheme)), FATAL)
          endif
         ! parse number of filling passes
           if (Control%Fill%fill_scheme(n) == EQUAL_BORROW) then
              if (parse(params,'hp', np) == 1) Control%Fill%npass_horiz(n) = np
              if (parse(params,'vp', np) == 1) Control%Fill%npass_vert (n) = np
           endif
      endif
   enddo

 ! print results
   if (mpp_pe() == mpp_root_pe()) then
           write (logunit,10) 'momentum', &
                 (trim(echo_advec_scheme(Control%scheme(m,-1))),m=1,2), Control%weight(-1)
           write (logunit,10) 'temperature', &
                 (trim(echo_advec_scheme(Control%scheme(m, 0))),m=1,2), Control%weight( 0)
      do n = 1, ntrace
           call get_tracer_names (MODEL_ATMOS, n, name)
           write (logunit,11) n, trim(name), &
                    ( trim(echo_advec_scheme(Control%scheme(m,n))),m=1,2), &
                                             Control%weight(n),            &
                                             Control%Fill%npass_horiz(n),  &
                                             Control%Fill%npass_vert(n)
        10 format (3x, a24, ', HORIZ=',a24, ', VERT=',a24, ', wt=',f7.3)
        11 format (i3, a24, ', HORIZ=',a24, ', VERT=',a24, ', wt=',f7.3, &
                   ', hp=',i2, ', vp=',i2)
      enddo
   endif

 ! check if fourth order centered horizontal scheme
   Control%do_mask4_vel = Control%scheme(1,-1) == FOURTH_CENTERED
   Control%do_mask4_tmp = .false.
   do n = 0, ntrace
     if (Control%scheme(1,n) == FOURTH_CENTERED) then
        Control%do_mask4_tmp = .true.
        exit
     endif
   enddo
 ! check if finite_volume horizontal scheme for temperature
     if (Control%scheme(1,0) == FINITE_VOLUME_LINEAR .or. &
         Control%scheme(1,0) == FINITE_VOLUME_PARABOLIC) then
         Control%do_finite_volume_tmp = .true.
     else
         Control%do_finite_volume_tmp = .false.
     endif
 ! check if finite_volume horizontal scheme for tracers
   Control%do_finite_volume_trs = .false.
   do n = 1, ntrace
     if (Control%scheme(1,n) == FINITE_VOLUME_LINEAR .or. &
         Control%scheme(1,n) == FINITE_VOLUME_PARABOLIC) then
         Control%do_finite_volume_trs = .true.
         exit
     endif
   enddo
 ! error checks
 ! many horizontal schemes are not support with this release
   do n = -1, ntrace
     if (Control%scheme(1,n) == FINITE_VOLUME_LINEAR    .or. &
         Control%scheme(1,n) == FINITE_VOLUME_PARABOLIC .or. &
         Control%scheme(1,n) == SECOND_CENTERED_WTS     .or. &
         Control%scheme(1,n) == FOURTH_CENTERED_WTS   ) then
         call error_mesg ('bgrid_advection_mod',  &
                   'advection scheme not supported', FATAL)
     endif
   enddo

! initialize code sections for performance timing 

  if (do_clock_init) then
    id_tmp  = mpp_clock_id ('BGRID: temperature advection', &
                        flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
    id_trs  = mpp_clock_id ('BGRID: tracer advection',      &
                        flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
    id_vel  = mpp_clock_id ('BGRID: momentum advection',    &
                        flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
    id_fill = mpp_clock_id ('BGRID: tracer borrowing',    &
                        flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
    do_clock_init = .false. 
  endif

!-----------------------------------------------------------------------

  end subroutine advection_init

!#######################################################################

 subroutine advection_end

  ! terminate vertical advection module
    call vert_advection_end

  ! terminate finite-volume advection module
  ! if (Control%do_finite_volume_tmp .or. Control%do_finite_volume_trs) then
  !     call horiz_advection_end
  ! endif   

  ! free up memory used (although not much)
    deallocate ( Control%scheme, Control%weight, &
                 Control%Fill%fill_scheme,       &       
                 Control%Fill%npass_horiz,       &       
                 Control%Fill%npass_vert         )       

 end subroutine advection_end

!#######################################################################
! converts character string names to integer parameters
! checks for all valid names, otherwise produces an error

 function set_advec_scheme ( scheme ) result ( advec_scheme )
 character(len=*), intent(in)  :: scheme
 integer                       :: advec_scheme

      if (uppercase(trim(scheme)) == 'SECOND_CENTERED') then
          advec_scheme = SECOND_CENTERED
      else if (uppercase(trim(scheme)) == 'FOURTH_CENTERED') then
          advec_scheme = FOURTH_CENTERED
      else if (uppercase(trim(scheme)) == 'FOURTH_CENTERED_WTS') then
          advec_scheme = FOURTH_CENTERED_WTS
      else if (uppercase(trim(scheme)) == 'SECOND_CENTERED_WTS') then
          advec_scheme = SECOND_CENTERED_WTS
      else if (uppercase(trim(scheme)) == 'FINITE_VOLUME_LINEAR') then
          advec_scheme = FINITE_VOLUME_LINEAR
      else if (uppercase(trim(scheme)) == 'FINITE_VOLUME_PARABOLIC') then
          advec_scheme = FINITE_VOLUME_PARABOLIC
      else if (uppercase(trim(scheme)) == 'NONE') then
          advec_scheme = NONE
      else
          call error_mesg ('bgrid_advection_mod',  &
               'invalid advection scheme, '//uppercase(trim(scheme)), FATAL)
      endif

 end function set_advec_scheme
!-----------------------------------------------
! converts integer parameter values to character string names
! checks for all valid schemes, otherwise produces an error

 function echo_advec_scheme ( scheme_number ) result ( scheme_name )
 integer, intent(in)  :: scheme_number
 character(len=24)    :: scheme_name

   select case (scheme_number)
       case(NONE)
            scheme_name = 'none                    '
       case(SECOND_CENTERED)
            scheme_name = 'second_centered         '
       case(FOURTH_CENTERED)
            scheme_name = 'fourth_centered         '
       case(SECOND_CENTERED_WTS)
            scheme_name = 'second_centered_wts     '
       case(FOURTH_CENTERED_WTS)
            scheme_name = 'fourth_centered_wts     '
       case(FINITE_VOLUME_LINEAR)
            scheme_name = 'finite_volume_linear    '
       case(FINITE_VOLUME_PARABOLIC)
            scheme_name = 'finite_volume_parabolic '
       case default
            call error_mesg  ('bgrid_advection_mod', &
                     'invalid advection scheme number', FATAL)
   end select

 end function echo_advec_scheme

!#######################################################################
!         routines to minimize negative tracer values 
!        conservative horizontal and vertical borrowing
!#######################################################################

subroutine horiz_borrow (Hgrid, dt, dpde, var, var_dt, mask, iters)

!-----------------------------------------------------------------------
!
!       removes negative values by borrowing from horizontal
!       neighbors.  (usually for specific humidity or tke )
!
!-----------------------------------------------------------------------
   type(horiz_grid_type), intent(inout)                      :: Hgrid
   real, intent(in)                                          :: dt
   real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpde
   real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: var
   real, intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: var_dt
   real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:), optional :: mask
integer, intent(in),    dimension(:),                       optional :: iters
!-----------------------------------------------------------------------
!  ( note:   0.0 < flt <= 0.125 )
   real, parameter :: flt  = 0.125
   real, parameter :: flt4 = flt*0.25
   real, parameter :: rmin = 0.0
!-----------------------------------------------------------------------

   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,size(var,3)) ::  &
                                hew, hns, radp, few, fns, rcur, rdif

   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: rew, rns
   real, dimension(Hgrid%jlb:Hgrid%jub) :: areax, areay

integer :: i, j, k, n, is, ie, js, je, nlev, ntrace, knt, mxknt(size(var,4))
!-----------------------------------------------------------------------

   mxknt = 1;  if (present(iters)) mxknt = iters
   if (maxval(mxknt) == 0) return

   nlev   = size(var,3)
   ntrace = size(var,4)

   is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
   js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

!-----------------------------------------------------------------------
!------ compute flux coeff common to all variables -------

   do j = js-1, je
     areax(j) = Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j)
     areay(j) = Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j+1)
   enddo

   do k = 1, nlev
   do j = js-1, je
   do i = is-1, ie
     hew(i,j,k) = areax(j)*(dpde(i,j,k)+dpde(i+1,j,k))
     hns(i,j,k) = areay(j)*(dpde(i,j,k)+dpde(i,j+1,k))
   enddo
   enddo
   enddo
 ! apply step-mountain mask?
   if (present(mask)) then
       do k = 1, nlev
       do j = js-1, je
       do i = is-1, ie
         hew(i,j,k) = hew(i,j,k)*mask(i,j,k)*mask(i+1,j,k)
         hns(i,j,k) = hns(i,j,k)*mask(i,j,k)*mask(i,j+1,k)
       enddo
       enddo
       enddo
   endif

   do k = 1, nlev
   do j = js, je
     radp(:,j,k) = flt4/(Hgrid%Tmp%area(j)*dpde(:,j,k))
   enddo
   enddo

!-----------------------------------------------------------------------
!---- tracer loop -----
!---- store current value of tracer in rcur ----

   do n = 1, ntrace

      rcur(:,:,:) = var(:,:,:,n) + var_dt(:,:,:,n)*dt

!---- iteration loop -----

      few = 0.0;  fns = 0.0

   do knt = 1, mxknt(n)

!--------------2-nd order lat/lon contributions-------------------------

!     --- do borrowing where adjacent values have opposite sign ---
!            but do not turn off fluxes previously turned on

      do k = 1, nlev
         do j = js-1, je
         do i = is-1, ie
            ! east-west fluxes
            if ((rcur(i+1,j,k) <  rmin .and. rcur(i,j,k) >= rmin) .or. &
                (rcur(i+1,j,k) >= rmin .and. rcur(i,j,k) <  rmin))     &
            few(i,j,k) = hew(i,j,k)
            rew(i,j) = (rcur(i+1,j,k)-rcur(i,j,k))*few(i,j,k)
            ! north-south fluxes
            if ((rcur(i,j+1,k) <  rmin .and. rcur(i,j,k) >= rmin) .or. &
                (rcur(i,j+1,k) >= rmin .and. rcur(i,j,k) <  rmin))     &
            fns(i,j,k) = hns(i,j,k)
            rns(i,j) = (rcur(i,j+1,k)-rcur(i,j,k))*fns(i,j,k)
         enddo
         enddo

         do j = js, je
         do i = is, ie
            rdif(i,j,k) = (rew(i,j)-rew(i-1,j)+ &
                           rns(i,j)-rns(i,j-1))*radp(i,j,k)
         enddo
         enddo
      enddo

     !---- halo update (assumed that halo size = 1) ----

      call update_halo (Hgrid, TEMP, rdif)

     !---- update current value of tracer ----

      rcur(:,:,:) = rcur(:,:,:) + rdif(:,:,:)

!-----------------------------------------------------------------------
   enddo ! iteration loop
!-----------------------------------------------------------------------

     !---- return the tendency -----

      var_dt(:,:,:,n) = (rcur(:,:,:) - var(:,:,:,n)) / dt

!-----------------------------------------------------------------------
   enddo ! tracer loop
!-----------------------------------------------------------------------

end subroutine horiz_borrow

!#######################################################################

   subroutine vert_borrow (dt, dpde, var, var_dt, iters)

!-----------------------------------------------------------------------
!
!  This removes negative specific humidity/mixing ratios by borrowing
!  from the grid boxes immediately above and below. If not enough
!  is available to fill the negative then a negative will remain.
!
!-----------------------------------------------------------------------
   real, intent(in)    :: dt, dpde(:,:,:), var(:,:,:,:)
   real, intent(inout) :: var_dt(:,:,:,:)
integer, intent(in), optional :: iters(:)
!-----------------------------------------------------------------------
   real, dimension(size(var,3)) :: deficit, surplus, rdpdt
   real :: divid, ratio_dn, ratio_up, var_dp
   integer :: i, j, k, n, m, nlev, num_iters(size(var,4)), num, num_var
!-----------------------------------------------------------------------

   num_iters = 4; if (present(iters)) num_iters = iters
   if (maxval(num_iters) == 0) return

   num_var = size(var,4)
   if (num_var == 0) return
   nlev = size(var,3)

!---- variable loop and iteration loop -----

 do n = 1, num_var
 do m = 1, num_iters(n)

!---- existing negatives will not be corrected ----

   do j = 1, size(var,2)
   do i = 1, size(var,1)

   do k = 1, nlev
      var_dp = (var(i,j,k,n)+dt*var_dt(i,j,k,n))*dpde(i,j,k)
      deficit(k) = min(var_dp,0.0)
      surplus(k) = max(var_dp,0.0)
      rdpdt(k) = 1./(dpde(i,j,k)*dt)
   enddo

!---- top level ----

   if (deficit(1) < 0.0) then
      divid = max(-surplus(2)/deficit(1),1.0)
      ratio_dn = surplus(2)/divid
      var_dt(i,j,1,n) = var_dt(i,j,1,n) + ratio_dn*rdpdt(1)
      var_dt(i,j,2,n) = var_dt(i,j,2,n) - ratio_dn*rdpdt(2)
      surplus(2) = max(surplus(2)-ratio_dn,0.)
   endif

!---- interior levels ----

   do k = 2, nlev-1
     if (deficit(k) < 0.0) then
       divid = max(-(surplus(k-1)+surplus(k+1))/deficit(k),1.0)
       ratio_up = surplus(k-1)/divid
       ratio_dn = surplus(k+1)/divid
       var_dt(i,j,k,n)   = var_dt(i,j,k,n)+(ratio_up+ratio_dn)*rdpdt(k)
       var_dt(i,j,k-1,n) = var_dt(i,j,k-1,n)-ratio_up*rdpdt(k-1)
       var_dt(i,j,k+1,n) = var_dt(i,j,k+1,n)-ratio_dn*rdpdt(k+1)
       surplus(k+1) = max(surplus(k+1)-ratio_dn,0.)
     endif
   enddo

!---- bottom level ----

   if (deficit(nlev) < 0.0) then
      divid = max(-surplus(nlev-1)/deficit(nlev),1.0)
      ratio_up = surplus(nlev-1)/divid
      var_dt(i,j,nlev  ,n) = var_dt(i,j,nlev  ,n)+ratio_up*rdpdt(nlev)
      var_dt(i,j,nlev-1,n) = var_dt(i,j,nlev-1,n)-ratio_up*rdpdt(nlev-1)
   endif

 enddo
 enddo

 enddo
 enddo

!-----------------------------------------------------------------------

   end subroutine vert_borrow

!#######################################################################

subroutine compute_etadot_vel ( Hgrid, Vgrid, Masks, few, fns, etadot )
!subroutine compute_etadot_vel ( Hgrid, Vgrid, Masks, res, few, fns, etadot )

type(horiz_grid_type), intent(in) :: Hgrid
type (vert_grid_type), intent(in) :: Vgrid
type (grid_mask_type), intent(in) :: Masks
!real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:)   :: res
real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns
real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: etadot

real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: sdiv
integer :: i, j, k, is, ie, js, je, nlev

! this code will not work with the eta coordinate
! the variable res need to be added as an argument
  if (.not.Masks%sigma) call error_mesg  ('bgrid_advection_mod', &
      'etadot cannot be recomputed for the eta coordinate', FATAL)


  is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
  js = Hgrid%Vel%js;  je = Hgrid%Vel%je
  nlev = size(few,3)

  etadot = 0.
  sdiv   = 0.

  do k = 1, nlev

      ! vertical integral of divergence at velocity points
      ! only need etadot for velocity compute domain
      do j = js, je
      do i = is, ie
         sdiv(i,j) = sdiv(i,j) + ((few(i+1,j,k)+fns(i,j+1,k))-  &  
                                  (few(i  ,j,k)+fns(i,j  ,k)))  &
                                 *Hgrid%Vel%rarea(j)
      enddo   
      enddo   

      if (k < nlev) etadot(:,:,k+1) = -sdiv(:,:)

   enddo

   ! finally include d(ps)/dt term
   do k = 2, nlev
     !etadot(:,:,k) = etadot(:,:,k) + sdiv(:,:)*res(:,:)*Vgrid%eta(k))
      etadot(:,:,k) = etadot(:,:,k) + sdiv(:,:)*Vgrid%eta(k)
      if (.not.Masks%sigma) etadot(:,:,k) = etadot(:,:,k)*Masks%Vel%mask(:,:,k)
   enddo

end subroutine compute_etadot_vel

!#######################################################################

end module bgrid_advection_mod

