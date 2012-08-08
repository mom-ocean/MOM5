! ============================================================================
! root uptake module
! ============================================================================
module uptake_mod

#include "../shared/debug.inc"

use constants_mod, only: PI
use fms_mod, only : write_version_number
use soil_tile_mod, only : &
     soil_tile_type, soil_pars_type, max_lev, psi_wilt
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN

public :: uptake_init

public :: darcy2d_uptake, darcy2d_uptake_solver
public :: darcy2d_uptake_lin, darcy2d_uptake_solver_lin
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'uptake',&
    version     = '$Id: uptake.F90,v 17.0 2009/07/21 03:03:04 fms Exp $',&
    tagname     = '$Name: siena_201207 $'

! values for internal soil uptake option selector
integer, parameter ::   &
     UPTAKE_LINEAR         = 1, &
     UPTAKE_DARCY2D        = 2, &
     UPTAKE_DARCY2D_LIN    = 3

! ==== module variables ======================================================
logical :: module_is_initialized =.FALSE.
integer :: num_l ! # of water layers
real    :: dz    (max_lev)    ! thicknesses of layers
real    :: zfull (max_lev)

contains

! ============================================================================
subroutine uptake_init(num_l_in, dz_in, zfull_in)
  integer, intent(in) :: num_l_in ! # of layers
  real   , intent(in) :: &
       dz_in(:), &  ! layer thickness
       zfull_in(:)  ! layer centers

  call write_version_number(version, tagname)
  module_is_initialized =.TRUE.

  num_l = num_l_in
  dz    = dz_in
  zfull = zfull_in
end subroutine uptake_init


! ============================================================================
! given soil and root parameters, calculate the flux of water toward root
! per unit root length, and its derivative w.r.t. xylem water potential
subroutine darcy2d_flow (psi_x, psi_soil, K_sat, psi_sat, b, K_r, r_r, R, eps, u, du, psi_root)
  real, intent(in) :: &
       psi_x,    & ! xylem water potential, m
       psi_soil, & ! soil water potential, m
       K_sat,    & ! saturated soil hydraulic conductivity, kg/(m2 s)
       psi_sat,  & ! saturates soil water potential, m
       b,        & ! power of soil moisture characteristic function
       K_r,      & ! root membrane permeability per unit area, kg/(m3 s)
       r_r,      & ! radius of root, m
       R,        & ! characteristic radial half-distance between roots, m
       eps         ! requested precision in terms of psi, m

  real, intent(out) :: &
       u,        & ! uptake, kg/(m s)
       du,       & ! derivative of uptake w.r.t psi_x, kg/(m2 s)
       psi_root    ! water potential at the root/soil interface, m

  ! ---- local constants
  integer :: max_iter = 50 ! max number of iterations
  ! ---- local vars
  real :: u_soil, u_root ! water flows through soil and root skin, respectively, kg/(m s)
  real :: f ! u_soil - u_root difference; we look for f(psi_root) = 0
  real :: df ! derivative of w.r.t psi_root
  real :: n
  real :: C_r ! 
  real :: K_s
  real :: K_root ! root membrane prmeability per unit length, kg/(m2 s)
  real :: pl, ph ! brackets of the solution
  real :: psi_root0 ! previous guess for root water potential
  real :: dpsi
  integer :: iter 

  C_r=2*PI/(log(R/r_r))
  n = -(1+3/b)
  K_root = 2*PI*r_r*K_r

  ! set up values bracketing the solution for psi_root, so that
  ! f(pl)<0 and f(ph)>0, where f = u_soil - u_root
  if(psi_soil>psi_x) then
     pl = psi_soil; ph = psi_x
  else
     ph = psi_soil; pl = psi_x
  endif
  ! choose initial values of f and df so that we always do bisection on the
  ! first iteration. 
  ! That means that our first approximation is (psi_soil+psi_x)/2 -- in future,
  ! modify to get something better
  f=1; df=0; psi_root=pl 

  do iter = 1, max_iter
     psi_root0=psi_root
     if (((psi_root-ph)*df-f)*((psi_root-pl)*df-f)>0) then
        ! Newton step would throws us out of range, do bisection
        psi_root = (pl+ph)/2
     else
        ! do Newton step
        psi_root = psi_root - f/df
     endif
     dpsi=psi_root-psi_root0

     ! calculate flux difference and its derivative
     u_soil = C_r*K_sat*&
          (psi_sat/n* &
          (  (min(psi_soil,psi_sat)/psi_sat)**n   &
            -(min(psi_root,psi_sat)/psi_sat)**n ) &
          + max(0.0, psi_soil - psi_sat)          &
          - max(0.0, psi_root - psi_sat)          )
     u_root = K_root*(psi_root-psi_x)
     
     f=u_soil-u_root
     df=-C_r*K_sat*(min(psi_root,psi_sat)/psi_sat)**(n-1)-K_root

     ! update brackets so that they still enclose the root
     if(f>0) then
        ph = psi_root
     else
        pl = psi_root
     endif

     if(abs(dpsi)<eps) exit
  enddo

  u = u_root; 
  ! calcilate derivalive of u w.r.t psi_x
  K_s = C_r*K_sat*(min(psi_root,psi_sat)/psi_sat)**(n-1)
  du = -K_root*K_s/(K_root+K_s)

end subroutine 


! ============================================================================
! given water potential of roots and soil, and array of vegetation factors in
! the uptake formula, calculate the total soil water uptake, its derivative
! w.r.t. root water potential, and optionally the vertical distribution of the
! uptake.
! NOTE that is we use one-way uptake option then U(psi_root) is continuous, but
! DUDpsi_root is not
subroutine darcy2d_uptake ( soil, psi_x0, VRL, K_r, r_r, uptake_oneway, &
     uptake_from_sat, uptake, duptake)
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: &
       psi_x0, &   ! water potential inside roots (in xylem) at zero depth, m
       VRL(:), &   ! Volumetric Root Length (root length per unit volume), m/m3
       K_r,    &   ! permeability of the root skin per unit area, kg/(m3 s)
       r_r         ! radius of the roots, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but 
                   ! never loose it to the soil
       uptake_from_sat   ! if false, uptake from saturated soil is prohibited
  real, intent(out) :: &
       uptake(:), & ! water uptake by roots
       duptake(:)   ! derivative of water uptake w.r.t. psi_root
  ! ---- local constants
  real, parameter :: eps = 1e-4 ! tolerance for psi
  ! ---- local vars
  integer :: l
  real :: psi_x     ! water potential inside roots (psi_x0+z), m
  real :: psi_soil  ! water potential of soil, m
  real :: psi_sat   ! saturation soil water potential, m
  real :: k_sat     ! hyraulic conductivity of saturated soil, kg/(m2 s)
  real :: R         ! characteristic half-distance between roots, m
  
  real :: u         ! water uptake by roots at the current layer, kg/(m2 s)
  real :: du        ! derivative of u w.r.t. root water potential
  real :: psi_r


  ! calculate some hydraulic properties common for all soil layers
  psi_sat = soil%pars%psi_sat_ref/soil%pars%alpha
  k_sat   = soil%pars%k_sat_ref*soil%pars%alpha**2

  if(is_watch_point())then
     write(*,*)'##### darcy2d_uptake input #####'
     __DEBUG3__(psi_x0,psi_sat,K_sat)
     __DEBUG2__(K_r,r_r)
  endif
  ! calculate soil water supply and its derivative
  uptake = 0; duptake = 0
  do l = 1, num_l
     psi_x    = psi_x0+zfull(l)
     psi_soil = soil%psi(l)
     if (VRL(l) > 0) then
        R     = 1.0/sqrt(PI*VRL(l)) ! characteristic half-distance between roots, m
     else
        R     = 1.0 ! the value doesn't matter since uptake is 0 anyway 
     endif

     if ( soil%prog(l)%ws > 0 ) &
          cycle ! skip layers with ice
     if ( uptake_oneway.and.psi_x > soil%psi(l) ) &
          cycle ! skip layers where roots would loose water
     if ( .not.(uptake_from_sat).and.psi_soil >= psi_sat ) &
          cycle ! skip layers where the soil is saturated

     ! calculates soil term of uptake expression
     call darcy2d_flow (psi_x, psi_soil, K_sat, psi_sat, soil%pars%chb, K_r, r_r, R, eps, u, du, psi_r)

     ! scale by volumetric root length and thickness of layer to get total uptake 
     ! from the current soil layer
     uptake(l)  = VRL(l)*dz(l)*u ; duptake(l) = VRL(l)*dz(l)*du
     if(is_watch_point()) then
        write(*,'(a,i2.2,100(2x,a,g))')'level=',l, &
             'VRL=', VRL(l), 'R=', R,&
             'psi_x=', psi_x, 'psi_r=', psi_r, 'psi_soil=', psi_soil, &
             'U=',u,&
             'z=', zfull(l)
     endif
  enddo
end subroutine darcy2d_uptake


! =============================================================================
! for Darcy-flow uptake, find the root water potential such to satisfy actual 
! uptake by the vegetation. 
subroutine darcy2d_uptake_solver (soil, vegn_uptk, VRL, K_r, r_r, uptake_oneway, &
     uptake_from_sat, uptake, n_iter)
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)  :: &
       vegn_uptk, & ! uptake requested by vegetation, kg/(m2 s)
       VRL(:),    & ! volumetric root length, m/m3
       K_r,       & ! root membrane permeability per unit area, kg/(m3 s)
       r_r          ! root radius, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but 
                    ! never loose it to the soil
       uptake_from_sat ! if false, uptake from saturated soil is prohibited
  real,    intent(out) :: uptake(:) ! soil water uptake, by layer
  integer, intent(out) :: n_iter ! # of iterations made, for diagnostics only

  real :: uptake_tot

  call uptake_solver_K(soil, vegn_uptk, VRL, K_r, r_r, uptake_oneway, &
     uptake_from_sat, uptake, n_iter, darcy2d_uptake)

  ! since the numerical solution is not exact, adjust the vertical profile 
  ! of uptake to ensure that the sum is equal to transpiration exactly
  uptake_tot = sum(uptake(:))
  uptake(:) = uptake(:)+(vegn_uptk-uptake_tot)/sum(dz(:))*dz(:) 
  
end subroutine darcy2d_uptake_solver

! =============================================================================
! kernel of the uptake solver: given the input and a subroutine that calculates 
! the uptake vertical profile for given water potential at the surface, returns
! a soulution
subroutine uptake_solver_K (soil, vegn_uptk, VRL, K_r, r_r, uptake_oneway, &
     uptake_from_sat, uptake, n_iter, uptake_subr)
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)  :: &
       vegn_uptk, & ! uptake requested by vegetation, kg/(m2 s)
       VRL(:),    & ! volumetric root length, m/m3
       K_r,       & ! root membrane permeability per unit area, kg/(m3 s)
       r_r          ! root radius, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but 
                    ! never loose it to the soil
       uptake_from_sat ! if false, uptake from saturated soil is prohibited
  real,    intent(out) :: uptake(:) ! vertical distribution of soil uptake
  integer, intent(out) :: n_iter ! # of iterations made, for diagnostics only

  interface 
     subroutine uptake_subr ( soil, psi_x0, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, uptake, duptake)
     use soil_tile_mod, only : soil_tile_type
       type(soil_tile_type), intent(in) :: soil
       real, intent(in) :: &
            psi_x0, &   ! water potential inside roots (in xylem) at zero depth, m
            VRL(:), &   ! Volumetric Root Length (root length per unit volume), m/m3
            K_r,    &   ! permeability of the root skin per unit area, kg/(m3 s)
            r_r         ! radius of the roots, m
       logical, intent(in) :: &
            uptake_oneway, & ! if true, then the roots can only take up water, but 
            ! never loose it to the soil
            uptake_from_sat   ! if false, uptake from saturated soil is prohibited
       real, intent(out) :: &
            uptake(:), & ! water uptake by roots
            duptake(:)   ! derivative of water uptake w.r.t. psi_root
     end subroutine uptake_subr
  end interface

  ! ---- local constants
  ! parameters of uptake equation solver:
  integer, parameter :: max_iter = 50    ! max number of iterations
  real   , parameter :: eps      = 1e-10 ! end condition

  ! ---- local vars
  real :: xl,xh,x2,f,DfDx, incr
  real :: duptake(size(uptake))
  integer :: i

  n_iter = 0

  xl = psi_wilt
  xh = maxval(soil%psi(1:num_l)-zfull(1:num_l))

  ! find the lower upper boundary of the interval that contains solution
  incr = 100.0 ! inital psi increment for the lower bracket search
  do i = 1,20
     call uptake_subr ( soil, xl, VRL, K_r, r_r, uptake_oneway, uptake_from_sat, &
          uptake, duptake )
     if (sum(uptake)>=vegn_uptk) exit
     xl = xl-incr; incr=incr*2
  enddo
  if (sum(uptake) <= vegn_uptk) then
     ! Uptake is still smaller than vegn_uptake (that is,
     ! transpiration). We got as close to actual solution as possible.
     if (is_watch_point()) then
        write(*,*)'###### failed to reach lower bracket of uptake #####'
        __DEBUG2__(xl,uptake)
     endif
     return
  endif
  
  ! find upper boundary of the interval that contains solution
  incr = 1.0 ! inital psi increment for the upper bracket search
  do i = 1,20
     call uptake_subr ( soil, xh, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, uptake, duptake)
     if (sum(uptake)<=vegn_uptk) exit
     xh = xh+incr; incr = incr*2
  enddo
  if (sum(uptake)>= vegn_uptk) then
     ! Could not reach the psi_root high enough for uptake from soil to be 
     ! smaller than the transpiration: this can happen when transpiration is 
     ! negative (due to numerics) and uptake is one-way
     return
  endif

  ! the root is bracketed, find it
  x2 = (xl+xh)/2
  call uptake_subr ( soil, x2, VRL, K_r, r_r, uptake_oneway, &
       uptake_from_sat, uptake, duptake )
  f = sum(uptake) - vegn_uptk
  DfDx = sum(duptake)
  do n_iter = 1, max_iter
     ! check if we already reached the desired precision
     if(abs(f)<eps)exit
           
     if (is_watch_point()) then
        write(*,*)'##### solution iteration iter=',n_iter
        __DEBUG5__(f,DfDx,xl,xh,x2)
        __DEBUG2__((x2-xl)*DfDx,(x2-xh)*DfDx)
     endif
           
     if (((x2-xl)*DfDx-f)*((x2-xh)*DfDx-f)>0) then
        ! the Newton-Raphson step would throw us out of the bonds of the interval,
        ! so we do the bisection
        x2 = (xl+xh)/2
        if(is_watch_point()) write(*,*) 'did bisection'
     else
        x2 = x2-f/DfDx
        if(is_watch_point()) write(*,*) 'did Newton-Raphson step'
     endif
           
     call uptake_subr ( soil, x2, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, uptake, duptake)
     f = sum(uptake) - vegn_uptk
     DfDx = sum(duptake)
     
     if (f>0) then
        xl = x2
     else
        xh = x2
     endif
     
     if(is_watch_point()) then
        write(*,*)'#### After iteration',n_iter
        __DEBUG2__(vegn_uptk,sum(uptake))
     endif
  enddo

end subroutine uptake_solver_K

! ============================================================================
! given soil and root parameters, calculate the flux of water toward root
! per unit root length, and its derivative w.r.t. xylem water potential
! this version calculates fluxes linearized around psi_root0
subroutine darcy2d_flow_lin (psi_x, psi_soil, psi_root0, K_sat, psi_sat, b, K_r, &
     r_r, R, u, du, psi_root)
  real, intent(in) :: &
       psi_x,    & ! xylem water potential, m
       psi_soil, & ! soil water potential, m
       psi_root0,& ! value of psi_root we linearize around, m
       K_sat,    & ! saturated soil hydraulic conductivity, kg/(m2 s)
       psi_sat,  & ! saturates soil water potential, m
       b,        & ! power of soil moisture characteristic function
       K_r,      & ! root membrane permeability per unit area, kg/(m3 s)
       r_r,      & ! radius of root, m
       R           ! characteristic radial half-distance between roots, m
  real, intent(out) :: &
       u,        & ! uptake, kg/(m s)
       du,       & ! derivative of uptake w.r.t psi_x, kg/(m2 s)
       psi_root    ! water potential at the root/soil interface, m

  ! ---- local vars
  real :: u_soil0 ! flux through soil for psi_root = psi_root0
  real :: du_soil ! its derivative w.r.t. psi_root
  real :: C_r !
  real :: n
  real :: K_root  ! root membrane prmeability per unit length, kg/(m2 s)

  C_r=2*PI/(log(R/r_r))
  n = -(1+3/b)
  K_root = 2*PI*r_r*K_r

  ! calculate flux through soil for psi_root = psi_root0
  u_soil0 = C_r*K_sat*&
       (psi_sat/n* &
          (  (min(psi_soil ,psi_sat)/psi_sat)**n   &
            -(min(psi_root0,psi_sat)/psi_sat)**n ) &
          + max(0.0, psi_soil  - psi_sat)          &
          - max(0.0, psi_root0 - psi_sat)          )
  ! and its derivative w.r.t. psi_root at psi_root0
  du_soil=-C_r*K_sat*(min(psi_root0,psi_sat)/psi_sat)**(n-1)

  ! flux through soil+membrane
  u  = K_root/(-du_soil+K_root)*(u_soil0+du_soil*(psi_x-psi_root0))
  ! and its derivative w.r.t. psi_x
  du = K_root/(-du_soil+K_root)*du_soil
  ! water potential at the root-soil interface
  psi_root = psi_x + u/K_root
end subroutine 


! ============================================================================
subroutine darcy2d_uptake_lin ( soil, psi_x0, VRL, K_r, r_r,uptake_oneway, &
    uptake_from_sat, u, du )
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: &
       psi_x0,    & ! water potential inside roots (in xylem) at zero depth, m
       VRL(:),    & ! Volumetric Root Length (root length per unit volume), m/m3
       K_r,       & ! permeability of the root membrane per unit area, kg/(m3 s)
       r_r          ! radius of fine roots, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but 
                   ! never loose it to the soil
       uptake_from_sat   ! if false, uptake from saturated soil is prohibited
  real, intent(out) :: &
       u(:), &      ! layer-by-layer distribution of uptake, kg/(m2 s)
       du(:)        ! derivative of u w.r.t. root water potential, kg/(m3 s)
  ! ---- local vars
  integer :: k
  real :: psi_x     ! water potential inside roots (psi_x0+z), m
  real :: psi_soil  ! water potential of soil, m
  real :: psi_sat   ! saturation soil water potential, m
  real :: K_sat     ! hyraulic conductivity of saturated soil, kg/(m2 s)
  real :: R         ! characteristic half-distance between roots, m
  
  real :: psi_root  ! water potential at the root/soil interface, m
  real :: psi_root0 ! initial guess of psi_root, m


  ! calculate some hydraulic properties common for all soil layers
  psi_sat = soil%pars%psi_sat_ref/soil%pars%alpha
  K_sat   = soil%pars%k_sat_ref*soil%pars%alpha**2

  u = 0; du = 0
  do k = 1, num_l
     psi_x    = psi_x0 + zfull(k)
     psi_soil = soil%psi(k)
     psi_root0= soil%psi(k) ! change it later to prev. time step value
     if (VRL(k)>0) then
        R     = 1.0/sqrt(PI*VRL(k)) ! characteristic half-distance between roots, m
     else
        R     = 1.0 ! the value doesn't matter since uptake is 0 anyway (no roots) 
     endif
     if ( soil%prog(k)%ws > 0 ) &
          cycle ! skip layers with ice
     if ( uptake_oneway.and.psi_x > soil%psi(k) ) &
          cycle ! skip layers where roots would loose water
     if ( .not.(uptake_from_sat).and.psi_soil >= psi_sat ) &
          cycle ! skip layers where the soil is saturated

     ! calculates soil term of uptake expression
     call darcy2d_flow_lin (psi_x, psi_soil, psi_root0, K_sat, psi_sat, soil%pars%chb, &
          K_r, r_r, R, u(k), du(k), psi_root)

     ! scale by volumetric root length and thickness of layer to get total 
     ! uptake from the current soil layer
     u(k)  = VRL(k)*dz(k)*u(k)
     du(k) = VRL(k)*dz(k)*du(k)
  enddo

end subroutine

! ============================================================================
subroutine darcy2d_uptake_solver_lin ( soil, vegn_uptk, VRL, K_r, r_r, &
     uptake_oneway, uptake_from_sat, uptake, n_iter )
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: &
       vegn_uptk, & ! uptake requested by vegetation, kg/(m2 s)
       VRL(:),    & ! Volumetric Root Length (root length per unit volume), m/m3
       K_r,       & ! permeability of the root membrane per unit area, kg/(m3 s)
       r_r          ! radius of fine roots, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but 
                   ! never loose it to the soil
       uptake_from_sat   ! if false, uptake from saturated soil is prohibited
  real, intent(out) :: &
       uptake(:)    ! layer-by-layer distribution of uptake, kg/(m2 s)
  integer, intent(out) :: n_iter ! # of iterations made, for diagnostics only

  ! ---- local vars
  integer :: k
  real :: psi_x0    ! water potential inside roots (in xylem) at zero depth, m
  
  real :: uptake_tot  ! total water uptake by roots, kg/(m2 s)
  real :: Duptake_tot ! derivative of water uptake w.r.t. psi_root, kg/(m3 s)
  real :: du(size(uptake)) ! derivative of u w.r.t. root water potential, kg/(m3 s)
  real :: d_psi_x     ! change of the xylem potential, m

  if (uptake_oneway) then
     call uptake_solver_K(soil, vegn_uptk, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, uptake, n_iter, darcy2d_uptake_lin )
     
     ! since the numerical solution is not exact, adjust the vertical profile 
     ! of uptake to ensure that the sum is equal to transpiration exactly
     uptake_tot = sum(uptake(:))
     uptake(:) = uptake(:)+(vegn_uptk-uptake_tot)/sum(dz(:))*dz(:) 
  else
     psi_x0 = 0.0
     call darcy2d_uptake_lin ( soil, psi_x0, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, uptake, du )
     uptake_tot = sum(uptake)
     Duptake_tot = sum(du) 
     
     if (duptake_tot/=0) then
        d_psi_x = (vegn_uptk - uptake_tot)/duptake_tot
        do k = 1,num_l
           uptake(k) = uptake(k) + d_psi_x*du(k)
        enddo
     else
        uptake(:) = 0
     endif
     n_iter = 0 ! because we didn't do any iterations
  endif

end subroutine


end module uptake_mod
