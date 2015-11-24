
module bgrid_horiz_diff_mod

!=======================================================================
!
!                    linear horizontal mixing
!               with option for any order accuracy
!
!=======================================================================

use mpp_mod, only: input_nml_file 
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_masks_mod      , only: grid_mask_type
use bgrid_prog_var_mod   , only: prog_var_type, var_init
use bgrid_halo_mod       , only: update_halo, vel_flux_boundary, &
                                 EAST, NORTH, NOPOLE, POLEONLY,  &
                                 TEMP, UWND, VWND
use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID

use         fms_mod, only:  error_mesg, FATAL, write_version_number,         &
                            open_namelist_file, check_nml_error, close_file, &
                            mpp_clock_id, mpp_clock_begin, mpp_clock_end,    &
                            MPP_CLOCK_SYNC, CLOCK_MODULE, file_exist,        &
                            mpp_pe, mpp_root_pe, uppercase, stdlog
use   constants_mod, only:  RADIUS

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_names, get_number_tracers

implicit none
private

!--------- public interfaces ----------

public   horiz_diff, horiz_diff_init

!-----------------------------------------------------------------------
!  namelist bgrid_horiz_diff_nml

!  damp_scheme_wind   Determines how horizontal damping coefficients
!  damp_scheme_temp    vary with latitude.  Possible values are: 0,1,2,3.
!                  [0 = no damping; 1 = coeff are uniform with latitude;
!                       2 = coeffs vary as 1/(dx**2+dy**2);
!                       3 = coeffs vary as 1/(dx**2);
!                       4 = same as scheme 1 but with scheme 2 applied
!                       poleward of reflat]
!                      Notes: Schemes 2-4 provide increased damping in
!                      higher latitudes. Temperature and tracers are 
!                      controlled using "damp_scheme_temp".

   integer :: damp_scheme_temp = 1
   integer :: damp_scheme_wind = 1

!  reflat       Latitude cutoff (in degrees) at which increased
!               high latitude damping is applied.  Equatorward of
!               this latitude uniform damping (scheme=1) is applied;
!               poleward of this latitude enhanced damping (scheme=2)
!               is applied.  This variable is only used when
!               damp_scheme_wind=4 or damp_scheme_temp=4.

   real :: reflat = 85.
 
!  damp_order_wind      The horizontal damping order for momentum,
!  damp_order_temp       temperature, and default order for all 
!  damp_order_tracer     prognostic tracers. Only even numbers are allowed.

   integer :: damp_order_wind   = 4   ! use: 0, 2, 4,....
   integer :: damp_order_temp   = 4
   integer :: damp_order_tracer = 4

!  damp_coeff_wind      The horizontal damping coefficients for
!  damp_coeff_temp        momentum, temperature, and default value for
!  damp_coeff_tracer      all prognostic tracers. Only positive values
!                         are allowed.

   real    :: damp_coeff_wind   = 0.35  
   real    :: damp_coeff_temp   = 0.35  
   real    :: damp_coeff_tracer = 0.35

!  slope_corr_wind      The topography slope correction for horizontal
!  slope_corr_temp        damping of momentum, temperature, and default
!  slope_corr_tracer      for all prognostic tracers. 

   real, dimension(4) :: slope_corr_wind   = (/0.,0.,0.,0./)
   real, dimension(4) :: slope_corr_temp   = (/0.,0.,0.,0./)
   real, dimension(4) :: slope_corr_tracer = (/0.,0.,0.,0./)

namelist /bgrid_horiz_diff_nml/ damp_scheme_temp, damp_scheme_wind,      &
                    damp_order_wind, damp_order_temp, damp_order_tracer, &
                    damp_coeff_wind, damp_coeff_temp, damp_coeff_tracer, &
                    slope_corr_wind, slope_corr_temp, slope_corr_tracer, &
                    reflat

!-----------------------------------------------------------------------
!--------- private data ----------

! control parameters for horizontal damping
! set via the namelist or field table
 type hdiff_control_type
     integer, pointer :: order(:) =>NULL()
     real   , pointer :: coeff(:) =>NULL()
     real   , pointer :: slope(:,:) =>NULL()
     logical, pointer :: do_slope_adj(:) =>NULL()
     logical       :: do_damping, do_slope_adj_temp
     integer       :: damping_scheme_wind, damping_scheme_temp
     real, dimension(:),   pointer :: areahx =>NULL(), &
                                      areahy =>NULL(), &
                                      areavx =>NULL(), &
                                      areavy => NULL()
     real, dimension(:,:), pointer :: wth => NULL(), &
                                      wtv =>NULL()
 end type hdiff_control_type

 type(hdiff_control_type), save :: Control

 integer  :: nlev ! number of model levels

 character(len=128) :: version='$Id: bgrid_horiz_diff.F90,v 19.0 2012/01/06 19:54:01 fms Exp $'
 character(len=128) :: tagname='$Name: tikal $'
 logical :: do_log = .true.

!  timing data
 integer :: id_total
 logical :: do_clock_init = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine horiz_diff ( Hgrid, Masks, nplev, dt, dpde, pres, Var, Var_dt )

!-----------------------------------------------------------------------
!
!   Hgrid  = horizontal grid constants
!   Masks  = grid masking constants for eta coordinate
!   nplev  = number of "pure" pressure levels at the top of the model
!   dt     = adjustment time step
!   dpde   = pressure weight for model layers
!   pres   = pressure at full model levels
!   Var    = prognostic variables at the last updated time level
!   Var_dt = tendency of prognostic variables since the last
!            updated time level
!
!-----------------------------------------------------------------------

type(horiz_grid_type), intent(inout)  :: Hgrid
type (grid_mask_type),    intent(in)  :: Masks

integer,                intent(in)    :: nplev
  real,                 intent(in)    :: dt
  real,                 intent(in)    :: dpde(Hgrid%ilb:,Hgrid%jlb:,:),&
                                         pres(Hgrid%ilb:,Hgrid%jlb:,:)
  type (prog_var_type), intent(in)    :: Var
  type (prog_var_type), intent(inout) :: Var_dt

!-----------------------------------------------------------------------

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,       &
                  size(dpde,3)) :: hkew3, hkns3, hkew, hkns, &
                                   hcew, hcns, dat, vdat, hdac

  real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, &
                     size(dpde,3),3) :: pterms

  real    :: dt_inv, hsign, xcoeff, ycoeff
  integer :: i, j, k, n, is, ie, js, je, ntp

!=======================================================================
!  --- should damping be done ??? ----

   if ( .not. Control%do_damping    ) return
   if ( Control%damping_scheme_wind+Control%damping_scheme_temp == 0 ) return
   call mpp_clock_begin (id_total)

!-----------------------------------------------------------------------
!       --- check the horizontal dimensions of the input array ---

      nlev = size(dpde,3)

    if (size(dpde,1) /= Hgrid%isize .or.               &
        size(dpde,2) /= Hgrid%jsize ) call error_mesg  &
       ('bgrid_horiz_diff', 'input array has the wrong dimensions.', FATAL)

!   ---- time step related values ----

    dt_inv  = 1./dt

!   ---- initialize flux weights ----

    hkew3 = 0.0;  hkns3 = 0.0

!-----------------------------------------------------------------------
!-------------setup temperature and tracer damping ---------------------
!-----------------------------------------------------------------------
   ntp = count( Control%order(1:Var_dt%ntrace) > 0 )
   if (Control%order(0) > 0 .or. ntp > 0) then

      ! compute grid indiecs
      is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
      js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

      do j = js, je
         hdac(:,j,:) = Hgrid%Tmp%rarea(j)
      enddo
      do k = nplev+1, nlev
      do j = js, je
         hdac(:,j,k) = hdac(:,j,k) / (2.0*dpde(:,j,k))
      enddo
      enddo

!     ----- mass-weighted fluxes -----
      do j = js-1, je
         hkew3(:,j,1:nplev) = Control%areahx(j)
         hkns3(:,j,1:nplev) = Control%areahy(j)
      enddo

      do k = nplev+1, nlev
      do j = js-1, je
      do i = is-1, ie
         hkew3(i,j,k) = Control%areahx(j) * (dpde(i,j,k)+dpde(i+1,j,k))
         hkns3(i,j,k) = Control%areahy(j) * (dpde(i,j,k)+dpde(i,j+1,k))
      enddo
      enddo
      enddo

!     ----- mask out fluxes below step-mountain ----

    if (.not.Masks%sigma) then
       do k = 1, nlev
       do j = js-1, je
       do i = is-1, ie
          hkew3(i,j,k) = hkew3(i,j,k)*Masks%Tmp%mask(i,j,k)*Masks%Tmp%mask(i+1,j,k)
          hkns3(i,j,k) = hkns3(i,j,k)*Masks%Tmp%mask(i,j,k)*Masks%Tmp%mask(i,j+1,k)
       enddo
       enddo
       enddo
          hdac(:,js:je,:) = hdac(:,js:je,:) * Masks%Tmp%mask(:,js:je,:)
    endif

!-----------------------------------------------------------------------
!----- slope adjustment ------------------------------------------------

      if ( Control % do_slope_adj_temp ) then
         call slope_correction_init ( Hgrid, Masks, nplev, pres, pterms )
      else
         hcew = 0.0
         hcns = 0.0
      endif

!-----------------------------------------------------------------------
!------------------temperature damping----------------------------------

      if (Control%order(0) > 0) then

         ! setting up flux weights - limit to 1/8
         do k = 1, nlev
         do j = Hgrid%jlb, Hgrid%jub
            xcoeff = min(0.125,Control%coeff(0)*Control%wth(j,1))
            ycoeff = min(0.125,Control%coeff(0)*Control%wth(j,2))
            hkew(:,j,k) = hkew3(:,j,k) * xcoeff
            hkns(:,j,k) = hkns3(:,j,k) * ycoeff
         enddo
         enddo

         dat  = Var%t + dt*Var_dt%t
!        --- slope adjustent ---
         if ( Control%do_slope_adj(0) ) &
         call slope_correction ( Hgrid, Masks, nplev, Control%slope(:,0), &
                                 pterms, dat, hcew, hcns )
         call diff_mass (Hgrid, dat, hcew, hcns, hkew, hkns, hdac,  &
                         Control%order(0))
         hsign = 1.; if (mod(Control%order(0),4) == 0) hsign = -1.
         Var_dt%t = Var_dt%t + hsign * dt_inv * dat

      endif

!-----------------------------------------------------------------------
!------------ tracer damping (prognostic tracers only) ---------------

      if (ntp > 0) then

         do n = 1, Var_dt%ntrace

            ! setting up flux weights - limit to 1/8
            do k = 1, nlev
            do j = Hgrid%jlb, Hgrid%jub
               xcoeff = min(0.125,Control%coeff(n)*Control%wth(j,1))
               ycoeff = min(0.125,Control%coeff(n)*Control%wth(j,2))
               hkew(:,j,k) = hkew3(:,j,k) * xcoeff
               hkns(:,j,k) = hkns3(:,j,k) * ycoeff
            enddo
            enddo

            if (Control%order(n) == 0) cycle
            dat  = Var%r(:,:,:,n) + dt*Var_dt%r(:,:,:,n)
           !--- slope adjustent ---
            if ( Control%do_slope_adj(n) ) &
            call slope_correction ( Hgrid, Masks, nplev, Control%slope(:,n), &
                                    pterms, dat, hcew, hcns )
            call diff_mass (Hgrid, dat, hcew, hcns, hkew, hkns, &
                            hdac, Control%order(n))
            hsign = 1.; if (mod(Control%order(n),4) == 0) hsign = -1.
            Var_dt%r(:,:,:,n) = Var_dt%r(:,:,:,n) + hsign * dt_inv * dat
         enddo

      endif

   endif
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------setup momentum damping -------------------------------

   if (Control%order(-1) > 0) then

      is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
      js = Hgrid%Vel%js;  je = Hgrid%Vel%je

      call change_grid (Hgrid, TEMP_GRID, WIND_GRID, &
                        dpde(:,:,nplev+1:nlev), vdat(:,:,nplev+1:nlev))
      call update_halo (Hgrid, UWND, vdat(:,:,nplev+1:nlev), &
                        halos=EAST+NORTH, flags=NOPOLE)

      do j = js, je
         hdac(:,j,:) = Hgrid%Vel%rarea(j)
      enddo
      do k = nplev+1, nlev
         hdac(:,js:je,k) = hdac(:,js:je,k) / (2.*vdat(:,js:je,k))
      enddo

!     ----- mass-weighted fluxes -----
      do j = js, je+1
         hkew3(:,j,1:nplev) = Control%areavx(j)
         hkns3(:,j,1:nplev) = Control%areavy(j)
      enddo

      do k = nplev+1, nlev
      do j = js, je+1
      do i = is, ie+1
         hkew3(i,j,k) = Control%areavx(j)*(vdat(i,j,k)+vdat(i-1,j,k))
         hkns3(i,j,k) = Control%areavy(j)*(vdat(i,j,k)+vdat(i,j-1,k))
      enddo
      enddo
      enddo

!     ----- mask out fluxes below step-mountain ----

    if (.not.Masks%sigma) then
       do k = 1, nlev
       do j = js, je+1
       do i = is, ie+1
          hkew3(i,j,k) = hkew3(i,j,k)*Masks%Vel%mask(i,j,k)*Masks%Vel%mask(i-1,j,k)
          hkns3(i,j,k) = hkns3(i,j,k)*Masks%Vel%mask(i,j,k)*Masks%Vel%mask(i,j-1,k)
       enddo
       enddo
       enddo
          hdac(:,js:je,:) = hdac(:,js:je,:) * Masks%Tmp%mask(:,js:je,:)
    endif

!-----------------------------------------------------------------------
!----- slope adjustment setup ------

    if ( Control % do_slope_adj(-1) ) then
!      ---- pressure at velocity points ----
       vdat(:,:,1:nplev) = pres(:,:,1:nplev)
       call change_grid ( Hgrid, TEMP_GRID, WIND_GRID, &
                          pres(:,:,nplev+1:nlev), vdat(:,:,nplev+1:nlev) )
       call update_halo (Hgrid, UWND, vdat(:,:,nplev+1:nlev), &
                                        halos=EAST+NORTH, flags=NOPOLE)

       call vel_slope_correction_init ( Hgrid, Masks, nplev, &
                                        Control%slope(:,-1), &
                                        vdat, pterms         )
    endif

!-----------------------------------------------------------------------

      ! setting up flux weights - limit to 1/8
      do k = 1, nlev
      do j = js, je+1
         xcoeff = min(0.125,Control%coeff(-1)*Control%wtv(j,1))
         ycoeff = min(0.125,Control%coeff(-1)*Control%wtv(j,2))
         hkew(:,j,k) = hkew3(:,j,k) * xcoeff
         hkns(:,j,k) = hkns3(:,j,k) * ycoeff
      enddo
      enddo
      ! zero-out cross polar fluxes
      call vel_flux_boundary (Hgrid, hkns)

!-----------------------------------------------------------------------
!-------------------------momentum damping------------------------------

       dat = Var%u + dt*Var_dt%u
      vdat = Var%v + dt*Var_dt%v
      if ( Control % do_slope_adj(-1) ) then
         if (Masks%sigma) then
             call diff_vel (Hgrid, dat,vdat, hkew,hkns, hdac,  &
                            Control%order(-1), pterms)
         else
             call diff_vel (Hgrid, dat,vdat, hkew,hkns, hdac,  &
                            Control%order(-1), pterms, Masks%Vel%kbot)
         endif
      else
         call diff_vel (Hgrid, dat,vdat, hkew,hkns, hdac,  &
                        Control%order(-1))
      endif
      hsign = 1.; if (mod(Control%order(-1),4) == 0) hsign = -1.
      Var_dt%u = Var_dt%u + hsign * dt_inv *  dat
      Var_dt%v = Var_dt%v + hsign * dt_inv * vdat

   endif

   call mpp_clock_end (id_total)

!-----------------------------------------------------------------------

 end subroutine horiz_diff

!#######################################################################
!#######################################################################

 subroutine diff_mass (Hgrid, rdat, hcew, hcns, hkew, hkns, hdac, &
                       order)

!----------------------------------------------------------------------
!
!        diff_mass is a private interface that performs multiple
!        2nd order lapacians.
!
!----------------------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)    :: order
   real   , intent(inout), dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: rdat
   real   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: hcew, hcns, &
                                              hkew, hkns, hdac

!----------------------------------------------------------------------
!  hcew, hcns are weighted corrections to the diffusive fluxes for
!  sloping sigma surfaces
!----------------------------------------------------------------------

   real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
              rew, rns
   integer :: i, j, k, n, is, ie, js, je, nordr

!-----------------------------------------------------------------------

   is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
   js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

!-----------------------------------------------------------------------
!------------loop for order of damping scheme-------------------------

   do n = 1, order/2

        do k = 1, nlev

!------------------contributions (fluxes) ------------------------------

     if ( n == 1 ) then
        do j = js-1, je
        do i = is-1, ie
           rew(i,j) = (rdat(i+1,j,k)-rdat(i,j,k)+hcew(i,j,k))*hkew(i,j,k)
           rns(i,j) = (rdat(i,j+1,k)-rdat(i,j,k)+hcns(i,j,k))*hkns(i,j,k)
        enddo
        enddo
     else
        do j = js-1, je
        do i = is-1, ie
           rew(i,j) = (rdat(i+1,j,k)-rdat(i,j,k))*hkew(i,j,k)
           rns(i,j) = (rdat(i,j+1,k)-rdat(i,j,k))*hkns(i,j,k)
        enddo
        enddo
     endif

!-----------------------------------------------------------------------

     do j = js, je
     do i = is, ie
        rdat(i,j,k)=(rew(i,j)-rew(i-1,j)+rns(i,j)-rns(i,j-1)) &
                    *hdac(i,j,k)
     enddo
     enddo

     enddo

!-----------------------------------------------------------------------
!---- update all halo rows ? ----
!   do not update on last pass, halos will updated in the main program

     if ( n < order/2 ) then
         call update_halo (Hgrid, TEMP, rdat)
     endif

!-----------------------------------------------------------------------

   enddo

!-----------------------------------------------------------------------

 end subroutine diff_mass

!#######################################################################

 subroutine diff_vel (Hgrid, udat, vdat, vkew, vkns, hdac, order, terms, kbot)

!----------------------------------------------------------------------
!
!        diff_vel is a private interface that performs multiple
!        2nd order lapacians for the momentum components.
!
!----------------------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)    :: order
   real   , intent(inout), dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: udat, vdat
   real   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: vkew, vkns, hdac
   real   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev, 3), optional :: terms
   integer, intent(in),    dimension(Hgrid%ilb:Hgrid%iub,  &
                                     Hgrid%jlb:Hgrid%jub), &
                                               optional :: kbot

!----------------------------------------------------------------------

   real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) ::  &
              uew, uns, vew, vns, dudp, dvdp, ucew, ucns, vcew, vcns
   integer :: i, j, k, n, is, ie, js, je, k1, k2

!-----------------------------------------------------------------------

   is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
   js = Hgrid%Vel%js;  je = Hgrid%Vel%je

!------------loop for order of damping scheme---------------------------

   do n = 1, order/2

      do k = 1, nlev

!------------------contributions (fluxes) ------------------------------

      if ( n == 1 .and. present(terms) ) then

        k1 = max(k-1,1)
        if (present(kbot)) then
           do j = Hgrid%jlb, Hgrid%jub
           do i = Hgrid%ilb, Hgrid%iub
              k2 = min(k+1,kbot(i,j))
              dudp(i,j) = (udat(i,j,k2)-udat(i,j,k1))*terms(i,j,k,1)
              dvdp(i,j) = (vdat(i,j,k2)-vdat(i,j,k1))*terms(i,j,k,1)
           enddo
           enddo
        else
           k2 = min(k+1,nlev)
           dudp(:,:) = (udat(:,:,k2)-udat(:,:,k1))*terms(:,:,k,1)
           dvdp(:,:) = (vdat(:,:,k2)-vdat(:,:,k1))*terms(:,:,k,1)
        endif

        do j = js, je+1
        do i = is, ie+1
!          ---- slope correction terms ----
           ucew(i,j) = (dudp(i,j)+dudp(i-1,j))*terms(i,j,k,2)
           ucns(i,j) = (dudp(i,j)+dudp(i,j-1))*terms(i,j,k,3)
           vcew(i,j) = (dvdp(i,j)+dvdp(i-1,j))*terms(i,j,k,2)
           vcns(i,j) = (dvdp(i,j)+dvdp(i,j-1))*terms(i,j,k,3)

           uew(i,j) = (udat(i,j,k)-udat(i-1,j  ,k)+ucew(i,j))*vkew(i,j,k)
           uns(i,j) = (udat(i,j,k)-udat(i  ,j-1,k)+ucns(i,j))*vkns(i,j,k)
           vew(i,j) = (vdat(i,j,k)-vdat(i-1,j  ,k)+vcew(i,j))*vkew(i,j,k)
           vns(i,j) = (vdat(i,j,k)-vdat(i  ,j-1,k)+vcns(i,j))*vkns(i,j,k)
        enddo
        enddo

      else

        do j = js, je+1
        do i = is, ie+1
           uew(i,j) = (udat(i,j,k)-udat(i-1,j  ,k))*vkew(i,j,k)
           uns(i,j) = (udat(i,j,k)-udat(i  ,j-1,k))*vkns(i,j,k)
           vew(i,j) = (vdat(i,j,k)-vdat(i-1,j  ,k))*vkew(i,j,k)
           vns(i,j) = (vdat(i,j,k)-vdat(i  ,j-1,k))*vkns(i,j,k)
        enddo
        enddo

      endif

!-----------------------------------------------------------------------
      do j = js, je
      do i = is, ie
        udat(i,j,k) = (uew(i+1,j  )-uew(i,j)+ &
                       uns(i  ,j+1)-uns(i,j))*hdac(i,j,k)
        vdat(i,j,k) = (vew(i+1,j  )-vew(i,j)+ &
                       vns(i  ,j+1)-vns(i,j))*hdac(i,j,k)
      enddo
      enddo

      enddo
!-----------------------------------------------------------------------
!---- update all halo rows ? ----
!   do not update on last pass, halos will updated in the main program

      if ( n < order/2 ) then
         call update_halo (Hgrid, UWND, udat)
         call update_halo (Hgrid, VWND, vdat)
      endif

!-----------------------------------------------------------------------

   enddo

!-----------------------------------------------------------------------

 end subroutine diff_vel

!#######################################################################

 subroutine horiz_diff_init ( Hgrid )

!-----------------------------------------------------------------------
!       initialization of horizontal damping coefficients
!-----------------------------------------------------------------------

 type(horiz_grid_type), intent(in) :: Hgrid ! horizontal grid constants

!-----------------------------------------------------------------------

real    :: eps = 1.e-6
integer :: j

 integer            :: n, nv, order, ntrace, unit, ierr, io, logunit
 real               :: coeff, slope(4)
 character(len=128) :: scheme, params, tname

!-----------------------------------------------------------------------
! read namelist
   if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=bgrid_horiz_diff_nml, iostat=io)
       ierr = check_nml_error(io,'bgrid_horiz_diff_nml')
#else
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read (unit, nml=bgrid_horiz_diff_nml, iostat=io, end=5)
          ierr = check_nml_error (io, 'bgrid_horiz_diff_nml')
       enddo
 5     call close_file (unit)
#endif
   endif

   logunit = stdlog()
   if (do_log) then
      call write_version_number(version, tagname)
      if (mpp_pe() == mpp_root_pe()) write (logunit, nml=bgrid_horiz_diff_nml)
      do_log = .false.
   endif

   call get_number_tracers ( MODEL_ATMOS, num_prog=ntrace )

   allocate ( Control%order        (-1:ntrace), &
              Control%coeff        (-1:ntrace), &
              Control%do_slope_adj (-1:ntrace), &
              Control%slope      (4,-1:ntrace)  )

 ! defaults
   Control%order = 4
   Control%coeff = .35
   Control%slope =  0.

 ! namelist arguments
   Control%order  (-1) = damp_order_wind
   Control%coeff  (-1) = damp_coeff_wind
   Control%slope(:,-1) = slope_corr_wind

   Control%order  (0) = damp_order_temp
   Control%coeff  (0) = damp_coeff_temp
   Control%slope(:,0) = slope_corr_temp

   Control%order  (1:ntrace) = damp_order_tracer
   Control%coeff  (1:ntrace) = damp_coeff_tracer
   Control%slope(:,1:ntrace) = spread(slope_corr_tracer,2,ntrace)

 ! set the damping scheme (check for errors below)
   Control % damping_scheme_wind = damp_scheme_wind
   Control % damping_scheme_temp = damp_scheme_temp

 ! process tracer table information for horizontal damping methods
   do n = 1, ntrace
      if (query_method('diff_horiz', MODEL_ATMOS, n, scheme, params)) then
          if (uppercase(trim(scheme)) /= 'LINEAR' .and. uppercase(trim(scheme)) /= 'NONE') &
              call error_mesg  ('bgrid_horiz_diff_mod',  &
                            'invalid damping method, '//uppercase(trim(scheme)), FATAL)
          if (parse(params,'order', order) == 1) Control%order(n) = order
          if (uppercase(trim(scheme)) == 'NONE') Control%order(n) = 0
          if (parse(params,'coeff', coeff) == 1) Control%coeff(n) = coeff
          nv = parse(params,'slope', slope)
          Control%slope (1:nv,n) = slope(1:nv)
      endif
      if (mpp_pe() == mpp_root_pe()) then
         call get_tracer_names (MODEL_ATMOS, n, tname)
         write (logunit,10) n, trim(tname), Control%order(n), Control%coeff(n), Control%slope(:,n)
      10 format (i3, a24, ', Order=',i2, ', Coeff=',f10.5, ', Slope=',4f10.5)
      endif
   enddo

 ! error checking
   do n = -1, ntrace
     !-- order
      if (Control%order(n) < 0 .or. mod(Control%order(n),2) /= 0) &
          call error_mesg ('bgrid_horiz_diff_mod', 'invalid damping order', FATAL)
     !-- non-dimension, normalized coefficient
      if (Control%coeff(n) < 0. .or. Control%coeff(n) > 1.) call error_mesg &
                   ('bgrid_horiz_diff_mod', 'invalid damping coeff', FATAL)
     !-- slope correction weights
      if (minval(Control%slope(:,n)) < 0. .or. maxval(Control%slope(:,n)) > 1.) &
            call error_mesg ('bgrid_horiz_diff_mod', 'invalid slope correction coeff', FATAL)
   enddo

 ! set flags
   Control%do_damping        = maxval(Control%order) > 0
   Control%do_slope_adj_temp = .false.
   do n = -1, ntrace
     Control%do_slope_adj(n) = maxval(Control%slope(:, n)) > 1.e-6
     ! set flag for temp OR tracer slope adjustment
     if (n >= 0 .and. Control%do_slope_adj(n)) Control%do_slope_adj_temp = .true.
   enddo

!-----------------------------------------------------------------------
!----- pre-compute metric terms ------

     allocate ( Control%areahx (Hgrid%jlb:Hgrid%jub), &
                Control%areahy (Hgrid%jlb:Hgrid%jub), &
                Control%areavx (Hgrid%jlb:Hgrid%jub), &
                Control%areavy (Hgrid%jlb:Hgrid%jub), &
                Control%wth    (Hgrid%jlb:Hgrid%jub,2), &
                Control%wtv    (Hgrid%jlb:Hgrid%jub,2)  )

!    ---- areas averaged along axes ----

     do j = Hgrid%Tmp%jsd, Hgrid%Tmp%jed-1
        Control%areahx(j) = 0.5*(Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j))
        Control%areahy(j) = 0.5*(Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j+1))
     enddo

     do j = Hgrid%Vel%jsd+1, Hgrid%Vel%jed
        Control%areavx(j) = 0.5*(Hgrid%Vel%area(j)+Hgrid%Vel%area(j))
        Control%areavy(j) = 0.5*(Hgrid%Vel%area(j)+Hgrid%Vel%area(j-1))
     enddo


!    ---- damping weight of x and y axis varies ----
!    ---- depending which damping scheme is used ----

     call damp_scheme_init ( Hgrid, WIND_GRID, &
                             Control%damping_scheme_wind, Control%wtv )
     call damp_scheme_init ( Hgrid, TEMP_GRID, &
                             Control%damping_scheme_temp, Control%wth )

! initialize code sections for performance timing 

  if (do_clock_init) then
     id_total = mpp_clock_id ('BGRID: horiz_diff (TOTAL)', &
                      flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE)
     do_clock_init = .false. 
  endif

!-----------------------------------------------------------------------

 end subroutine horiz_diff_init

!#######################################################################

 subroutine damp_scheme_init ( Hgrid, grid, scheme, weights )
 type(horiz_grid_type), intent(in) :: Hgrid
 integer, intent(in)  :: grid, scheme
 real,    intent(out) :: weights(Hgrid%jlb:Hgrid%jub,2)
 
 real, dimension(Hgrid%jlb:Hgrid%jub) :: dx2
 real    :: dxdy2eq, dx2eq, dy2, factor2
 integer :: j
 

!    ---- damping weight of x and y axis varies ----
!    ---- depending which damping scheme is used ----
!
!      scheme 1:   equal/constant
!      scheme 2:   function of diagonal grid distance
!      scheme 3:   function of x-axis grid distance
!

  select case (scheme)

!   ---- uniform damping ----
    case (1)
       weights = 0.125

    case (2:5)
       dx2=0.0
       if (grid == TEMP_GRID) then
       do j = Hgrid%jlb, Hgrid%jub-1
             dx2(j) = Hgrid%Tmp%dx(j)**2
       enddo
             dy2 = Hgrid%Tmp%dy**2
       else if (grid == WIND_GRID) then
       do j = Hgrid%jlb+1, Hgrid%jub
             dx2(j) = Hgrid%Vel%dx(j)**2
       enddo
             dy2 = Hgrid%Vel%dy**2
       else
          ! error condition needed
       endif

       ! function of diagonal distance (poleward of reflat)
       if (scheme == 2 .or. scheme == 4) then
           factor2 = 1./cos(reflat*acos(0.)/90.)**2
           dxdy2eq = RADIUS**2*(Hgrid%dlm**2+Hgrid%dph**2)
           weights(:,1) = 0.125*max(1.,dxdy2eq/(factor2*dx2+dy2))
           if (scheme == 2) then
               weights(:,2) = weights(:,1)
           else
               weights(:,2) = 0.125
           endif
       endif

       ! function of x-distance 
       if (scheme == 3 .or. scheme == 5) then
           dx2eq = (RADIUS*Hgrid%dlm*cos(reflat*acos(0.)/90.))**2
           where (dx2 /= 0.0) weights(:,1) = 0.125*max(1.,dx2eq/dx2)
           if (scheme == 3) then
               weights(:,2) = weights(:,1)
           else
               weights(:,2) = 0.125
           endif
       endif


    case default

       call error_mesg ('bgrid_horiz_diff_mod', 'invalid damping scheme', &
                         FATAL)

  end select

 end subroutine damp_scheme_init

!#######################################################################

 subroutine slope_correction_init ( Hgrid, Masks, nplev, pres, terms )

  type(horiz_grid_type), intent(in)  :: Hgrid
  type (grid_mask_type), intent(in)  :: Masks
  integer,               intent(in)  :: nplev
  real,                  intent(in)  :: pres (Hgrid%ilb:,Hgrid%jlb:,:)
  real,                  intent(out) :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  integer :: i, j, k, k1, k2, nlev

!  initialization of pressure terms for the sigma slope correction
!  these pressure terms do not change between mass variables
!  may want make weight a function of variable and/or level
!      USE ONE-HALF OF SPECIFIED WEIGHT AT LOWEST LEVEL

   nlev = size(pres,3)

   do k = nplev+1, nlev
     !--- reciprocal of vert gradient ---
      k1 = max(k-1,1)
      if (Masks%sigma) then
         k2 = min(k+1,nlev)
         terms(:,:,k,1) = 1.0/(pres(:,:,k2)-pres(:,:,k1))
      else
         do j = Hgrid%jlb, Hgrid%jub
         do i = Hgrid%ilb, Hgrid%iub
            k2 = min(k+1,Masks%Tmp%kbot(i,j))
            if (k2 > k1) then
              terms(i,j,k,1) = 1.0/(pres(i,j,k2)-pres(i,j,k1))
            else
              terms(i,j,k,1) = 0.0
            endif
         enddo
         enddo
      endif

!     --- horiz gradients (flip sign) ---
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub-1
           terms(i,j,k,2) = (pres(i,j,k)-pres(i+1,j,k))
           terms(i,j,k,3) = (pres(i,j,k)-pres(i,j+1,k))
      enddo
      enddo
   enddo

 end subroutine slope_correction_init

!#######################################################################

 subroutine vel_slope_correction_init ( Hgrid, Masks, nplev, weights, pres, terms )

  type(horiz_grid_type), intent(in)  :: Hgrid
  type (grid_mask_type), intent(in)  :: Masks
  integer,               intent(in)  :: nplev
  real,                  intent(in)  :: weights(4)
  real,                  intent(in)  :: pres (Hgrid%ilb:,Hgrid%jlb:,:)
  real,                  intent(out) :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  real    :: wt2
  integer :: i, j, k, k1, k2, ks, nlev, isd, ied, jsd, jed
  real :: dp(size(pres,1),size(pres,2))

!  initialization of pressure terms for the sigma slope correction
!  these pressure terms do not change between mass variables
!  may want make weight a function of variable and/or level
!      USE ONE-HALF OF SPECIFIED WEIGHT AT LOWEST LEVEL

   nlev = size(pres,3)
   isd = Hgrid%Vel%isd; ied = Hgrid%Vel%ied
   jsd = Hgrid%Vel%jsd; jed = Hgrid%Vel%jed

   do k = 1, nplev
      terms(:,:,k,:) = 0.0
   enddo

   do k = nplev+1, nlev
      k1 = max(k-1,1)
      if (Masks%sigma) then
          k2 = min(k+1,nlev)
          dp = (pres(:,:,k2)-pres(:,:,k1))
      else
          do j = Hgrid%jlb, Hgrid%jub
          do i = Hgrid%ilb, Hgrid%iub
             k2 = min(k+1,Masks%Vel%kbot(i,j))
             dp(i,j) = (pres(i,j,k2)-pres(i,j,k1))
          enddo
          enddo
      endif
     !--- reciprocal of vert gradient ---
      where (dp > 0.)
        terms(:,:,k,1) = 1.0/dp
      elsewhere
        terms(:,:,k,1) = 1.e30 ! these values should not be used where it counts
      endwhere

!     --- horiz gradients (flip sign) ---
      ks = max(1,k-nlev+4)
      wt2  = 0.5*weights(ks)
      do j = Hgrid%jlb+1, Hgrid%jub
      do i = Hgrid%ilb+1, Hgrid%iub
           terms(i,j,k,2) = wt2*(pres(i-1,j,k)-pres(i,j,k))
           terms(i,j,k,3) = wt2*(pres(i,j-1,k)-pres(i,j,k))
      enddo
      enddo
   enddo

 end subroutine vel_slope_correction_init

!#######################################################################

 subroutine slope_correction ( Hgrid, Masks, nplev, weights, terms, temp, cew, cns )

  type(horiz_grid_type), intent(in)  :: Hgrid
  type (grid_mask_type), intent(in)  :: Masks
  integer,               intent(in)  :: nplev
  real,                  intent(in)  :: weights(4)
  real,                  intent(in)  :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  real,                  intent(in)  :: temp (Hgrid%ilb:,Hgrid%jlb:,:)
  real,                  intent(out) :: cew  (Hgrid%ilb:,Hgrid%jlb:,:),&
                                        cns  (Hgrid%ilb:,Hgrid%jlb:,:)
  integer :: i, j, k, k1, k2, ks, nlev
  real :: wt2
  real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: dtdp

!  computes weighted-corrections for the slope of sigma surfaces
!  to east-west and north-south fluxes of field temp

 ! check for no correction
   do k = 1, 4
      if (weights(k) > 1.e-6) go to 10
   enddo
   cew = 0.0; cns = 0.0
   return

10 nlev = size(temp,3)

   do k = 1, nplev
      cew(:,:,k) = 0.0
      cns(:,:,k) = 0.0
   enddo

   do k = nplev+1, nlev
      k1 = max(k-1,1)
      if (Masks%sigma) then
         k2 = min(k+1,nlev) ! one-sided derivative at surface
         dtdp(:,:) = (temp(:,:,k2)-temp(:,:,k1))*terms(:,:,k,1)
      else
         do j = Hgrid%jlb, Hgrid%jub
         do i = Hgrid%ilb, Hgrid%iub
            k2 = min(k+1,Masks%Tmp%kbot(i,j))
            dtdp(i,j) = (temp(i,j,k2)-temp(i,j,k1))*terms(i,j,k,1)
         enddo
         enddo
      endif

      ks = max(1,k-nlev+4)
      wt2  = 0.5*weights(ks)
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub-1
           ! note: terms are grouped for reproducibility with previous version
           cew(i,j,k) = (dtdp(i,j)+dtdp(i+1,j))*(terms(i,j,k,2)*wt2)
           cns(i,j,k) = (dtdp(i,j)+dtdp(i,j+1))*(terms(i,j,k,3)*wt2)
      enddo
      enddo
   enddo

 end subroutine slope_correction

!#######################################################################

end module bgrid_horiz_diff_mod

