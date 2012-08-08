
module bgrid_horiz_adjust_mod

!-----------------------------------------------------------------------
!   This modules has interfaces for computing various 
!         horizontal adjustment processes.
!
!   1) momentum adjustment using coriolis and pressure gradient
!   2) mass fluxes (used for divergence and advection)
!   3) mass divergence
!   4) thermodynamic term (horizontal part of omega-alpha term)
!   4) pressure gradient (several options)
!   5) grad(p)/p term (used for energy conservation)
!   6) divergence damping
!
!-----------------------------------------------------------------------

use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, compute_pres_depth, &
                                 compute_geop_height, compute_pres_half
use bgrid_masks_mod      , only: grid_mask_type
use bgrid_halo_mod       , only: update_halo, TEMP
use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID, &
                                              UFLX_GRID, VFLX_GRID

use         constants_mod, only: OMEGA, RADIUS, RDGAS
use               fms_mod, only: error_mesg, FATAL

implicit none
private

public :: horiz_adjust_vel, horiz_adjust_mass, press_grad, &
          div_damping, compute_grad_pres, press_grad_fv

!-----------------------------------------------------------------------

   real, parameter :: RADIUS_INV = 1.0/RADIUS

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine horiz_adjust_vel ( Hgrid, Masks, dt, &
                              pgfew, pgfns, um, vm, udt, vdt,  &
                              alpha_implicit )

!-----------------------------------------------------------------------
!  update momentum tendencies using coriolis and pressure gradient
!
!    IN: Hgrid  = horizontal grid constants
!        Masks  = grid masking constants
!        dt     = time step
!        pgfew,
!         pgfns = pressure gradient force components (m/s2)
!        um,vm  = zonal and meridional wind components
!
! INOUT: udt,vdt = tendency of zonal and meridional wind components
!
! IN (opt): alpha_implicit = coefficient for coriolis/pgf time diff
!                               0.0 = explicit (not recommended)
!                               0.5 = trapezoidal implicit (default)
!                               1.0 = fully implicit
!
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in)      :: Hgrid
type (grid_mask_type), intent(in)      :: Masks
   real, intent(in)                    :: dt
   real, intent(in),    dimension(Hgrid%ilb:, Hgrid%jlb:, :) :: &
                                           pgfew, pgfns, um, vm
   real, intent(inout), dimension(Hgrid%ilb:, Hgrid%jlb:, :) :: &
                                                        udt, vdt
real, optional, intent(in)              :: alpha_implicit
!-----------------------------------------------------------------------

  real, dimension (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
      pgfu, pgfv, f0, fa, fb, cu, cv, up, vp, un, vn

    integer :: i, j, k
    real    :: alpha
!-----------------------------------------------------------------------

     alpha = 0.5; if (present(alpha_implicit)) alpha = alpha_implicit
     alpha = min(max(0.0,alpha),1.0)

!-----------------------------------------------------------------------
!-----------------update u and v (coriolis & pgf)-----------------------

      do k =            1, size(um,3)
      do j = Hgrid%Vel%js, Hgrid%Vel%je
      do i = Hgrid%Vel%is, Hgrid%Vel%ie

!     --------- pressure gradient force components---------------

         pgfu(i,j) = dt*pgfew(i,j,k)
         pgfv(i,j) = dt*pgfns(i,j,k)

!     ------------coriolis & curvature terms----------------------

         f0(i,j) = (um(i,j,k)*RADIUS_INV*Hgrid%tanphv(i,j)+  &
                              2.*OMEGA*Hgrid%sinphv(i,j))*dt
         fa(i,j) = f0(i,j)*(1.0-alpha)
         fb(i,j) = f0(i,j)*alpha

         cu(i,j) =  vm(i,j,k)*fa(i,j)
         cv(i,j) = -um(i,j,k)*fa(i,j)

!     ------------compute new u and v (coriolis & pgf)------------

         up(i,j) = pgfu(i,j) + cu(i,j) + um(i,j,k)
         vp(i,j) = pgfv(i,j) + cv(i,j) + vm(i,j,k)

         un(i,j) = ((fb(i,j)*vp(i,j)+up(i,j))/   &
                    (fb(i,j)*fb(i,j)+1.0))*Masks%Vel%mask(i,j,k)
         vn(i,j) = (vp(i,j)-fb(i,j)*un(i,j))*Masks%Vel%mask(i,j,k)

!     ---- return unfiltered tendencies with halos not updated ----

         udt(i,j,k) = udt(i,j,k) + (un(i,j)-um(i,j,k))/dt
         vdt(i,j,k) = vdt(i,j,k) + (vn(i,j)-vm(i,j,k))/dt
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine horiz_adjust_vel

!#######################################################################

 subroutine horiz_adjust_mass ( nplev, Hgrid, Masks,  &
                                u, v, dpde, cew, cns,        &
                                flew, flns, div, omgalf      )

!-----------------------------------------------------------------------
!  Computes:  mass fluxes, divergence, thermodynamic term (horiz part)
!
!    IN: nplev  = vertical index of uppermost pure pressure level
!                 (use nplev=0 for sigma models)
!        Hgrid  = horizontal grid constants
!        Masks  = grid masking constants
!        u, v   = prognostic variables for zonal and meridional wind
!        dpde   = pressure thickness of model layers
!        cew,
!          cns  = zonal and meridional components of grad(p)/p (no units)
!
! INOUT: flew,
!          flns = zonal, meridional mass fluxes (summation) (Pa-m2/s)
!
!   OUT: div    = mass divergence (Pa/s)
!        omgalf = horizontal part of omega-alpha term (Pa/s)
!
!-----------------------------------------------------------------------
integer, intent(in)                   :: nplev
type(horiz_grid_type), intent(inout)  :: Hgrid
type (grid_mask_type), intent(in)     :: Masks
  real, intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u, v, &
                                                        dpde, cew, cns
  real, intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: flew, flns
  real, intent(out),   dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: div, omgalf
!-----------------------------------------------------------------------

    real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::   &
                                  few, fns, udy, vdx, tew, tns, adpdxy

    integer :: i, j, k, is, ie, js, je


    is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
    js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      few = 0.0;  fns = 0.0

   do k = 1, size(u,3)

     !---- average mass weights for mass fluxes ----
      call change_grid ( Hgrid, TEMP_GRID, WIND_GRID, dpde(:,:,k), adpdxy )

!---- compute mass fluxes, divergence & horizontal omega-alpha term ----
!-------add in compute flux corrections for grid separation ------------
!---- sum input/output mass fluxes -------

      do j = js-1, je+1
         udy(:,j)=Hgrid%Vel%dy   *u(:,j,k)*adpdxy(:,j)   *Masks%Vel%mask(:,j,k)
         vdx(:,j)=Hgrid%Vel%dx(j)*v(:,j,k)*adpdxy(:,j)   *Masks%Vel%mask(:,j,k)
      enddo

      do j = js,   je
      do i = is-1, ie
          few(i,j)   = (udy(i,j)+udy(i,j-1))*0.5
          tew(i,j)   = few(i,j)*cew(i,j,k)
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
          fns(i,j)   = (vdx(i,j)+vdx(i-1,j))*0.5
          tns(i,j)   = fns(i,j)*cns(i,j,k)
      enddo
      enddo

!-----------------------------------------------------------------------
! ------ sum fluxes (output needed for advection) ------
! ------ (halos will be updated in advection) ------

      flew(:,:,k) = flew(:,:,k) + few
      flns(:,:,k) = flns(:,:,k) + fns

!---------------------------divergence----------------------------------

      do j = js, je
      do i = is, ie
         div(i,j,k) = ((few(i  ,j)+fns(i,j  ))-  &
                       (few(i-1,j)+fns(i,j-1)))  &
                        *Hgrid%Tmp%rarea(j)*Masks%Tmp%mask(i,j,k)
      enddo
      enddo

!-------------- horizontal part of omega-alpha -------------------------
!        ------ do not do for pure pressure levels -----

      if (k > nplev) then
         do j = js, je
         do i = is, ie
            omgalf(i,j,k)=(tew(i,j)+tew(i-1,j)+tns(i,j)+tns(i,j-1)) &
                          *0.50*Hgrid%Tmp%rarea(j)*Masks%Tmp%mask(i,j,k)
         enddo
         enddo
      else
            omgalf(:,:,k) = 0.0
      endif

!-----------------------------------------------------------------------
!---- end level loop -----

   enddo

!-----------------------------------------------------------------------

 end subroutine horiz_adjust_mass

!#######################################################################
!        Simmons and Burridge (1981) pressure gradient

subroutine press_grad ( Hgrid, Vgrid, Masks, fssl, tq,         &
                        dpde, wta, wtb, cew, cns, pgfew, pgfns )

!-----------------------------------------------------------------------
!    IN: Hgrid  = horizontal constants
!        Vgrid  = vertical constants
!        Masks  = grid masking constants
!        fssl   = geopotential height (m2/s2) at eta=1.
!        tq     = virtual temperature
!        dpde   = pressure thickness of model layers
!        wta,
!          wtb  = weights for computing geopotential height
!                 (same as weight for computing full pressures)
!        cew,
!          cns  = zonal and meridional components of grad(p)/p
!
!   OUT: pgfew,
!         pgfns = pressure gradient force components (m/s2)
!
!-----------------------------------------------------------------------
type(horiz_grid_type), intent (inout) :: Hgrid
type (vert_grid_type), intent (in)  :: Vgrid
type (grid_mask_type), intent (in)  :: Masks
   real, intent (in) , dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub) :: fssl
   real, intent (in),  dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       tq, dpde, wta, wtb, cew, cns
   real, intent (out), dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       pgfew, pgfns

!-----------------------------------------------------------------------
  real, dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                                       pew, pns, tew, tns

  real, dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,  &
                   size(tq,3)) :: fim, ppcew, ppcns

    integer :: i, j, k, is, ie, js, je
!-----------------------------------------------------------------------
!------------- integration of geopotential height-----------------------
!        fim = geopotential height at model levels

   if (Vgrid%nlev /= size(tq,3)) call error_mesg ('geop_height',  &
                              'wrong number of vertical levels', FATAL)

   if (Vgrid%nlev > 1) then

       call compute_geop_height (Vgrid, fssl, tq, wta, wtb, fim, &
                                 mask=Masks%Tmp%mask)

   else

       fim(:,:,:) = dpde(:,:,:) ! special case for shallow water model

   endif

!----------- lat/lon contributions to pressure gradient ----------------
!   -------- two loops: pressure only, pressure/sigma -------
!   ppcew,ppcns = zonal and meridional auxillary pressure gradient
!                 force components (i.e., not on velocity grid)

   is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
   js = Hgrid%Vel%js;  je = Hgrid%Vel%je

   ! compute at pure pressure levels
   do k =  1, Vgrid%nplev
     do j = js, je+1
     do i = is, ie
       ppcew(i,j,k) = fim(i+1,j,k)-fim(i,j,k)
     enddo
     enddo
     do j = js, je
     do i = is, ie+1
       ppcns(i,j,k) = fim(i,j+1,k)-fim(i,j,k)
     enddo
     enddo
   enddo

   ! compute at sigma/pressure levels
   do k = Vgrid%nplev+1, Vgrid%nlev
     do j = js, je+1
     do i = is, ie
         pew(i,j)   = fim(i+1,j,k)-fim(i,j,k)
         tew(i,j)   = RDGAS*0.50*(tq(i+1,j,k)+tq(i,j,k))
       ppcew(i,j,k) = pew(i,j)+tew(i,j)*cew(i,j,k)
     enddo
     enddo
     do j = js, je
     do i = is, ie+1
         pns(i,j)   = fim(i,j+1,k)-fim(i,j,k)
         tns(i,j)   = RDGAS*0.50*(tq(i,j+1,k)+tq(i,j,k))
       ppcns(i,j,k) = pns(i,j)+tns(i,j)*cns(i,j,k)
     enddo
     enddo
   enddo

!--------------compute pressure gradient force components---------------

   do k =  1, Vgrid%nlev
   do j = js, je
   do i = is, ie
      pgfew(i,j,k)=-0.50*(ppcew(i,j,k)+ppcew(i,j+1,k))*Hgrid%Vel%rdx(j)
      pgfns(i,j,k)=-0.50*(ppcns(i,j,k)+ppcns(i+1,j,k))*Hgrid%Vel%rdy
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------

 end subroutine press_grad

!#######################################################################
!      code for Lin (1997) finite volume pressure gradient

 subroutine press_grad_fv ( Hgrid, Vgrid, Masks, fssl, tq, phalf, &
                            wta, wtb, pgfew, pgfns                )       

!-----------------------------------------------------------------------
!    IN: Hgrid  = horizontal constants
!        Vgrid  = horizontal constants
!        Masks  = grid masking constants
!        fssl   = geopotential height (m2/s2) at eta=1.
!        tq     = virtual temperature
!        phalf  = pressure at model layer interfaces
!        wta,    
!          wtb  = weights for computing geopotential height
!                 (same as weight for computing full pressures)
!   OUT: pgfew,
!         pgfns = pressure gradient force components (m/s2)
!-----------------------------------------------------------------------
type(horiz_grid_type), intent (inout) :: Hgrid
type (vert_grid_type), intent (in)  :: Vgrid
type (grid_mask_type), intent (in)  :: Masks
   real, intent (in) , dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub) :: fssl 
   real, intent (in),  dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev+1) :: &
                                       phalf
   real, intent (in),  dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       tq, wta, wtb
   real, intent (out), dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       pgfew, pgfns

!-----------------------------------------------------------------------
  real, dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                                       few, pew, fns, pns
  real, dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,  &
                   size(tq,3)+1) :: fi, pi
  real   :: dp1, dp2
  integer :: i, j, k, is, ie, js, je

!------------- integration of geopotential height-----------------------
!        fi = geopotential height at half levels
!        pi = ln pressure at half levels

   if (Vgrid%nlev /= size(tq,3)) call error_mesg ('geop_height',  &
                              'wrong number of vertical levels', FATAL)

 ! initialize
     is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
     js = Hgrid%Vel%js;  je = Hgrid%Vel%je

 ! define ln pressure
     pi(:,:,2:Vgrid%nlev+1) = log(phalf(:,:,2:Vgrid%nlev+1))
     if (Vgrid%pzero) then
 ! extrapolate when pressure at top = 0
       ! following simmons & burridge scheme
       pi(:,:,1) = 2.*pi(:,:,2)*phalf(:,:,2)/(phalf(:,:,2)-phalf(:,:,1)) - &
                      pi(:,:,2) - 2.
     else
       pi(:,:,1) = log(phalf(:,:,1))
     endif

 ! compute geop height at model layer interfaces
 ! integrate up from surface
     fi(:,:,Vgrid%nlev+1) = fssl
  if (Masks%sigma) then
     do k = Vgrid%nlev, 1, -1
        fi(:,:,k) = fi(:,:,k+1) + RDGAS*tq(:,:,k)*(pi(:,:,k+1)-pi(:,:,k))
     enddo
  else
     do k = Vgrid%nlev, 1, -1
        where (Masks%Tmp%mask(:,:,k) < 0.5)
            fi(:,:,k) = Vgrid%fhalf(k)
        elsewhere
            fi(:,:,k) = fi(:,:,k+1) + RDGAS*tq(:,:,k)*(pi(:,:,k+1)-pi(:,:,k))
        endwhere
     enddo
  endif

 ! loop over levels
   do k = 1, Vgrid%nlev

   ! compute pressure gradient at intermediate points
   ! E/W component
     do j = js, je+1
     do i = is, ie
        dp1 = pi(i+1,j,k+1) - pi(i  ,j,k)
        dp2 = pi(i  ,j,k+1) - pi(i+1,j,k)
        few(i,j) = (dp1*(fi(i,j,k+1)-fi(i+1,j,k)) + dp2*(fi(i,j,k)-fi(i+1,j,k+1)))
        pew(i,j) = dp1+dp2
     enddo
     enddo
   ! N/S component
     do j = js, je
     do i = is, ie+1
        dp1 = pi(i,j+1,k+1) - pi(i,j  ,k)
        dp2 = pi(i,j  ,k+1) - pi(i,j+1,k)
        fns(i,j) = (dp1*(fi(i,j,k+1)-fi(i,j+1,k)) + dp2*(fi(i,j,k)-fi(i,j+1,k+1)))
        pns(i,j) = dp1+dp2
     enddo
     enddo

   ! compute pressure gradient at velocity points
     do j = js, je
     do i = is, ie
        pgfew(i,j,k) = (few(i,j)+few(i,j+1))/(pew(i,j)+pew(i,j+1))*Hgrid%Vel%rdx(j)
        pgfns(i,j,k) = (fns(i,j)+fns(i+1,j))/(pns(i,j)+pns(i+1,j))*Hgrid%Vel%rdy
     enddo
     enddo

   enddo

 end subroutine press_grad_fv

!#######################################################################
!   compute grad P term consistent with Simmons and Burridge

 subroutine compute_grad_pres (Hgrid, nplev, phalf, dpde, &
                               wta, wtb, cew, cns)

!-----------------------------------------------------------------------
!   IN:  Hgrid  = horizontal constants
!        nplev  = vertical index of uppermost pure pressure level
!                 (use nplev=0 for sigma models)
!        phalf  = pressure at model layer interfaces
!        dpde   = pressure thickness of model layers
!        wta,
!          wtb  = weights for computing geopotential height
!                 (same as weight for computing full pressures)
!  OUT:  cew,
!          cns  = zonal and meridional components of grad(p)/p (no units)
!-----------------------------------------------------------------------

   type(horiz_grid_type), intent(in)  :: Hgrid
   integer,               intent(in)  :: nplev
   real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) ::  &
                                       phalf, dpde, wta, wtb
   real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: cew, cns

   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                                   pewa, pewb, pnsa, pnsb
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                                        size(dpde,3)) :: adpdx, adpdy
   integer :: i, j, k

 ! special case for one level model
   if (size(dpde,3) == 1) then
       cew(:,:,:) = 0.0
       cns(:,:,:) = 0.0
       return
   endif

 ! term vanishes on pressure levels 
   do k = 1, nplev
       cew(:,:,k) = 0.0
       cns(:,:,k) = 0.0
   enddo

!---- compute pressure depth at uflx and vflx points ----

   call change_grid (Hgrid, TEMP_GRID, UFLX_GRID, dpde(:,:,:), adpdx)
   call change_grid (Hgrid, TEMP_GRID, VFLX_GRID, dpde(:,:,:), adpdy)


! --- compute pressure differences ----

   do j = Hgrid%Tmp%js,   Hgrid%Tmp%je+1
   do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
      pewa(i,j) = phalf(i+1,j,nplev+1)-phalf(i,j,nplev+1)
   enddo
   enddo

   do j = Hgrid%Tmp%js-1, Hgrid%Tmp%je
   do i = Hgrid%Tmp%is,   Hgrid%Tmp%ie+1
      pnsa(i,j) = phalf(i,j+1,nplev+1)-phalf(i,j,nplev+1)
   enddo
   enddo

! --- compute grad pressure term at model levels ----

   do k = nplev+1, size(phalf,3)-1

      !---- east-west contribution ----
       do j = Hgrid%Tmp%js,   Hgrid%Tmp%je+1
       do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
           pewb(i,j) = phalf(i+1,j,k+1)-phalf(i,j,k+1)
           cew(i,j,k) = (wta(i+1,j,k)+wta(i,j,k)) * pewa(i,j) + &
                        (wtb(i+1,j,k)+wtb(i,j,k)) * pewb(i,j)
           cew(i,j,k) = 0.5 * cew(i,j,k) / adpdx(i,j,k)
           pewa(i,j) = pewb(i,j)
       enddo
       enddo
      !---- north-south contribution ----
       do j = Hgrid%Tmp%js-1, Hgrid%Tmp%je
       do i = Hgrid%Tmp%is,   Hgrid%Tmp%ie+1
           pnsb(i,j) = phalf(i,j+1,k+1)-phalf(i,j,k+1)
           cns(i,j,k) = (wta(i,j+1,k)+wta(i,j,k)) * pnsa(i,j) + &
                        (wtb(i,j+1,k)+wtb(i,j,k)) * pnsb(i,j)
           cns(i,j,k) = 0.5 * cns(i,j,k) / adpdy(i,j,k)
           pnsa(i,j) = pnsb(i,j)
       enddo
       enddo
   enddo

 end subroutine compute_grad_pres

!#######################################################################

 subroutine div_damping (Hgrid, Vgrid, Masks, dt, coeff, dpde, div,  &
                         u_dt, v_dt)

!-----------------------------------------------------------------------
!
!   Computes divergence damping tendency for momentum fields
!
!    IN: Hgrid  = horizontal constants
!        Vgrid  = horizontal constants
!        Masks  = grid masking constants
!        dt     = time step
!        coeff  = coefficient for divergence damping
!        dpde   = model layer pressure thickness
!        div    = mass divergence (Pa/s)
!
! INOUT: u_dt,
!          v_dt = zonal and meridional momentum tendencies (m/s2)
!
!-----------------------------------------------------------------------

type(horiz_grid_type), intent (inout)  :: Hgrid
type (vert_grid_type), intent (in)  :: Vgrid
type (grid_mask_type), intent (in)  :: Masks
 real, intent (in)                  :: dt, coeff
 real, intent (in) ,   dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dpde, div
 real, intent (inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u_dt, v_dt

!-----------------------------------------------------------------------

  real, dimension (Hgrid%jlb:Hgrid%jub) :: scoeff
  real, dimension (Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub) :: dew, dns, adpdxy
  real, dimension (Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub,size(div,3)) :: adiv

  real    :: dcoeff
  integer :: i, j, k, is, ie, hs, he, vs, ve, ks, ke
  logical :: fourth_order
!-----------------------------------------------------------------------

  is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
  vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je
  hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
  ks = 1
  ke = size(div,3)
  fourth_order = .false.
  if (coeff < 0.0) fourth_order = .true.
  scoeff = abs(coeff)
  ! maximum damping at the sub-pole row
  ! if (vs == Hgrid%Vel%jsg) scoeff(vs) = 1.0
  ! if (ve == Hgrid%Vel%jeg) scoeff(ve) = 1.0

  if (fourth_order) then
     ! additional 5-pt smoother when fourth-order
     scoeff = -0.50*(scoeff/8.)**2/dt
     do k = ks, ke
        do j = hs-1, he
        do i = is-1, ie
           dew(i,j) = div(i+1,j,k)-div(i,j,k)
           dns(i,j) = div(i,j+1,k)-div(i,j,k)
        enddo
        enddo

        do j = hs, he
        do i = is, ie
           adiv(i,j,k) = (dew(i,j)+dns(i,j)-dew(i-1,j)-dns(i,j-1))
        enddo
        enddo
     enddo

     call update_halo ( Hgrid, TEMP, adiv )
  else
     ! second order
     scoeff = 0.50*(scoeff/8.)/dt
     do k = ks, ke
     do j = hs-1, he+1
        adiv(:,j,k) = div(:,j,k)
     enddo
     enddo
  endif

! second-order scheme (only the gradient)
  do k = ks, ke
     do j = vs, ve+1
     do i = is, ie
        dew(i,j) = adiv(i+1,j,k)-adiv(i,j,k)
     enddo
     enddo
     do j = vs, ve
     do i = is, ie+1
        dns(i,j) = adiv(i,j+1,k)-adiv(i,j,k)
     enddo
     enddo

     call change_grid ( Hgrid, TEMP_GRID, WIND_GRID, dpde(:,:,k), adpdxy )

     do j = vs, ve
     do i = is, ie
        dcoeff = scoeff(j) / adpdxy(i,j) * Masks%Vel%mask(i,j,k)
        u_dt(i,j,k) = u_dt(i,j,k) + (dew(i,j)+dew(i,j+1)) * dcoeff * Hgrid%Vel%dx(j)
        v_dt(i,j,k) = v_dt(i,j,k) + (dns(i,j)+dns(i+1,j)) * dcoeff * Hgrid%Vel%dx(j)
     enddo
     enddo
  enddo

!-----------------------------------------------------------------------

 end subroutine div_damping

!#######################################################################

end module bgrid_horiz_adjust_mod

