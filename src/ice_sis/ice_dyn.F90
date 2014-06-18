!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from         !
! Hunke and Dukowicz (JPO 1990, H&D hereafter) -Mike Winton (Michael.Winton) !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_dyn_mod

  use mpp_domains_mod, only: mpp_update_domains, BGRID_NE
  use constants_mod,   only: grav, pi
  use ice_grid_mod,    only: Domain, isc, iec, im, jsc, jec, isd, ied, jsd, jed, jm
  use ice_grid_mod,    only: dtw, dte, dts, dtn, dxt, dxv, dyt, dyv, cor, wett, wetv
  use ice_grid_mod,    only: t_on_uv, t_to_uv, dTdx, dTdy, dt_evp, evp_sub_steps
  use ice_grid_mod,    only: dydx, dxdy
  use ice_grid_mod,    only: reproduce_siena_201303
  use ice_thm_mod,     only: DI, DS, DW

  implicit none
  private

  public :: ice_dyn_param, ice_dynamics, strain_angle, ice_strength, sigI, sigII
  public :: blturn

  logical         :: SLAB_ICE = .false.  ! should we do old style GFDL slab ice?
  !
  ! parameters for calculating water drag and internal ice stresses
  !
  real            :: p0 = 2.75e4         ! pressure constant (Pa)
  real            :: c0 = 20.0           ! another pressure constant
  real            :: cdw = 3.24e-3       ! ice/water drag coef.
  real            :: blturn = 25.0       ! air/water surf. turning angle (NH) 25
  real, parameter :: EC = 2.0            ! yield curve axis ratio
  real, parameter :: EC2I = 1.0/(EC*EC)
  real, parameter :: MIV_MIN =  1.0      ! min ice mass to do dynamics (kg/m^2)

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_dyn_param - set ice dynamic parameters                                   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_dyn_param(p0_in, c0_in, cdw_in, wd_turn_in, slab_ice_in)
    real,    intent(in)   :: p0_in, c0_in, cdw_in, wd_turn_in
    logical, intent(in)   :: slab_ice_in

    p0 = p0_in
    c0 = c0_in
    cdw = cdw_in
    blturn = wd_turn_in
    SLAB_ICE = slab_ice_in
  end subroutine ice_dyn_param

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! set_strn - calculate generalized orthogonal coordinate strain tensor         !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine set_strn(ui, vi, strn11, strn22, strn12) ! ??? may change to do loop
    real, intent(in ), dimension(isd:ied,jsd:jed) :: ui, vi
    real, intent(out), dimension(isc:iec,jsc:jec) :: strn11, strn22, strn12
    real, allocatable, dimension(:,:), save :: fac1, fac2, fac3, fac4
    logical, save :: initialized = .false.
    integer       :: i, j

    if (.not. initialized) then
       allocate( fac1(isc:iec,jsc:jec), fac2(isc:iec,jsc:jec), &
                 fac3(isc:iec,jsc:jec), fac4(isc:iec,jsc:jec)  )
       do j = jsc, jec
          do i = isc, iec
             fac1(i,j) = (dtn(i,j)-dts(i,j))/dyt(i,j)
             fac2(i,j) = (dte(i,j)-dtw(i,j))/dxt(i,j)
             fac3(i,j) = 0.5*dyt(i,j)/dxt(i,j)
             fac4(i,j) = 0.5*dxt(i,j)/dyt(i,j)
          enddo
       enddo
       initialized = .true.
    end if

    do j = jsc, jec
       do i = isc, iec
          strn11(i,j) = (0.5*(ui(i,j)-ui(i-1,j)+ui(i,j-1)-ui(i-1,j-1))        &
                      + 0.25*(vi(i,j)+vi(i,j-1)+vi(i-1,j)+vi(i-1,j-1))*fac1(i,j))/dxt(i,j)
          strn22(i,j) = (0.5*(vi(i,j)-vi(i,j-1)+vi(i-1,j)-vi(i-1,j-1))        &
                      + 0.25*(ui(i,j)+ui(i,j-1)+ui(i-1,j)+ui(i-1,j-1))*fac2(i,j))/dyt(i,j)
          strn12(i,j) = fac3(i,j)*(0.5*(vi(i,j)/dyv(i,j)-vi(i-1,j)/dyv(i-1,j) &
                      + vi(i,j-1)/dyv(i,j-1)-vi(i-1,j-1)/dyv(i-1,j-1)))       &
                      + fac4(i,j)*(0.5*(ui(i,j)/dxv(i,j)-ui(i,j-1)/dxv(i,j-1) &
                      + ui(i-1,j)/dxv(i-1,j)-ui(i-1,j-1)/dxv(i-1,j-1)))
       enddo
    enddo

    return
  end subroutine set_strn

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_strength - magnitude of force on ice in plastic deformation              !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function ice_strength(hi, ci) ! ??? may change to do loop
!!$    real, dimension(isc:iec,jsc:jec), intent(in) :: hi, ci
    real, dimension(isc:,jsc:), intent(in) :: hi, ci
    real, dimension(isc:iec,jsc:jec)             :: ice_strength

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          ice_strength(i,j) = p0*hi(i,j)*ci(i,j)*exp(-c0*(1-ci(i,j)))
       enddo
    enddo

    return
  end function ice_strength

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_dynamics - take a single dynamics timestep with EVP subcycles            !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_dynamics(ci, hs, hi, ui, vi, sig11, sig22, sig12, uo, vo,       &
       fxat, fyat, sea_lev, fxoc, fyoc, fxic, fyic, fxco, fyco)
!!$    real, intent(in   ), dimension(isd:ied,jsd:jed) :: ci, hs, hi  ! ice properties
    real, intent(in   ), dimension(isd:,jsd:) :: ci, hs, hi  ! ice properties
    real, intent(inout), dimension(isd:ied,jsd:jed) :: ui, vi      ! ice velocity
    real, intent(inout), dimension(isd:ied,jsd:jed) :: sig11, sig22, sig12       ! stress tensor
    real, intent(in   ), dimension(isd:ied,jsd:jed) :: uo, vo      ! ocean velocity
!!$    real, intent(in   ), dimension(isc:iec,jsc:jec) :: fxat, fyat  ! air stress on ice
    real, intent(in   ), dimension(isc:,jsc:) :: fxat, fyat  ! air stress on ice
    real, intent(in   ), dimension(isd:ied,jsd:jed) :: sea_lev     ! sea level
    real, intent(  out), dimension(isc:iec,jsc:jec) :: fxoc, fyoc  ! ice stress on ocean
    real, intent(  out), dimension(isc:iec,jsc:jec) :: fxic, fyic  ! ice int. stress
    real, intent(  out), dimension(isc:iec,jsc:jec) :: fxco, fyco  ! coriolis force

    real, dimension(isc:iec,jsc:jec)    :: prs                    ! ice pressure
    real                                :: zeta, eta              ! bulk/shear viscosities
    real, dimension(isc:iec,jsc:jec)    :: strn11, strn12, strn22 ! strain tensor

    real,    dimension(isc:iec,jsc:jec) :: miv                 ! mass on v-points
    real,    dimension(isc:iec,jsc:jec) :: civ                 ! conc. on v-points
    complex                             :: rr                  ! linear drag coefficient
    real                                :: fxic_now, fyic_now  ! ice internal stress

    ! temporaries for ice stress calculation
    real                             :: del2, a, b, tmp
    real, dimension(isc:iec,jsc:jec) :: edt, mp4z, t0, t1, t2
    real, dimension(isc:iec,jsc:jec) :: f11, f22
    real, dimension(isc:iec,jsc:jec) :: sldx, sldy
    real, dimension(isc:iec,jsc:jec) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

    ! for velocity calculation
    real,    dimension(isc:iec,jsc:jec) :: dtmiv, rpart, fpart, uvfac
    complex                             :: newuv
    integer                             :: i,j,l

    fxoc = 0.0; fyoc = 0.0 ! zero these for summing later
    fxic = 0.0; fyic = 0.0
    fxco = 0.0; fyco = 0.0

    if (SLAB_ICE) then
       ui = uo; vi = vo;
       fxoc = fxat; fyoc = fyat;
       return
    end if

    if (evp_sub_steps==0) return;
    !
    ! sea level slope force
    !
    do j = jsc, jec
       do i = isc, iec
          sldx(i,j) = -dt_evp*grav*(0.5*(sea_lev(i+1,j+1)-sea_lev(i,j+1) &
                    + sea_lev(i+1,j)-sea_lev(i,j)))/dxv(i,j)
          sldy(i,j) = -dt_evp*grav*(0.5*(sea_lev(i+1,j+1)-sea_lev(i+1,j) &
                    + sea_lev(i,j+1)-sea_lev(i,j)))/dyv(i,j)
       enddo
    enddo

    ! put ice/snow mass and concentration on v-grid
    call t_to_uv(ci*(hi*DI+hs*DS), miv)
    call t_to_uv(ci, civ)
    !
    ! precompute prs, elastic timestep parameter, and linear drag coefficient
    !
    prs = ice_strength(hi(isc:iec,jsc:jec), ci(isc:iec,jsc:jec) )

    do j = jsc, jec
       do i = isc, iec
          if(dxt(i,j) < dyt(i,j) ) then
             edt(i,j) = (DI*dxt(i,j)*dxt(i,j)*ci(i,j)*hi(i,j))/(2*dt_evp)
          else
             edt(i,j) = (DI*dyt(i,j)*dyt(i,j)*ci(i,j)*hi(i,j))/(2*dt_evp)
          endif
       enddo
    enddo

    do j = jsc, jec
       do i = isc, iec
          if((wetv(i,j)>0.5) .and. miv(i,j) > MIV_MIN ) then ! values for velocity calculation (on v-grid)
             dtmiv(i,j) = dt_evp/miv(i,j)
          else
             ui(i,j) = 0.0
             vi(i,j) = 0.0
          endif
       enddo
    enddo

    do l=1,evp_sub_steps
       !
       ! calculate strain tensor for viscosities and forcing elastic eqn.
       if(reproduce_siena_201303) then
          call mpp_update_domains(ui, vi, Domain)
       else
          call mpp_update_domains(ui, vi, Domain, gridtype=BGRID_NE)
       endif
!rab       call mpp_update_domains(vi, Domain)
       !
       call set_strn(ui, vi, strn11, strn22, strn12)
       !
       ! calculate viscosities - how often should we do this ?
       !
       if (l>=1) then
          do j = jsc, jec
             do i = isc, iec
                del2 = (strn11(i,j)*strn11(i,j)+strn22(i,j)*strn22(i,j))*(1+EC2I)     &
                       +4*EC2I*strn12(i,j)*strn12(i,j) +2*strn11(i,j)*strn22(i,j)*(1-EC2I)  ! H&D eqn 9
                if(del2 > 4e-18 ) then
                   zeta = 0.5*prs(i,j)/sqrt(del2)
                else
                   zeta = 2.5e8*prs(i,j)
                endif

                if(zeta<4e8) zeta = 4e8 ! Hibler uses to prevent nonlinear instability

                eta = zeta*EC2I
                !
                ! some helpful temporaries
                !
                if(hi(i,j) > 0.0 ) then
                   mp4z(i,j) = -prs(i,j)/(4*zeta)
                   t0(i,j)   = 2*eta/(2*eta+edt(i,j))
                   tmp       = 1/(4*eta*zeta)
                   a         = 1/edt(i,j)+(zeta+eta)*tmp
                   b         = (zeta-eta)*tmp
                   t1(i,j)   = b/a
                   t2(i,j)   = a-b*b/a
                endif
             enddo
          enddo
       end if
       !
       ! timestep stress tensor (H&D eqn 21)
       !

       do j = jsc, jec
          do i = isc, iec
             if( (wett(i,j)>0.5) .and.(ci(i,j)*(DI*hi(i,j)+DS*hs(i,j))>MIV_MIN) ) then
                f11(i,j)   = mp4z(i,j)+sig11(i,j)/edt(i,j)+strn11(i,j)
                f22(i,j)   = mp4z(i,j)+sig22(i,j)/edt(i,j)+strn22(i,j)
                sig11(i,j) = (t1(i,j)*f22(i,j)+f11(i,j))/t2(i,j)
                sig22(i,j) = (t1(i,j)*f11(i,j)+f22(i,j))/t2(i,j)
                sig12(i,j) = t0(i,j)*(sig12(i,j)+edt(i,j)*strn12(i,j))
             else
                sig11(i,j) = 0.0
                sig22(i,j) = 0.0
                sig12(i,j) = 0.0 ! eliminate internal ice forces 
             endif
          enddo
       enddo

       call mpp_update_domains(sig11, Domain, complete=.false.)
       call mpp_update_domains(sig22, Domain, complete=.false.)
       call mpp_update_domains(sig12, Domain, complete=.true.)

       tmp1 = dTdy(sig12*dxt)
       tmp2 = dTdx(sig11*dyt)
       tmp3 = t_on_uv(sig12)
       tmp4 = t_on_uv(sig22)
       tmp5 = t_on_uv(sig11)
       tmp6 = dTdx(sig12*dyt)
       tmp7 = dTdy(sig22*dxt)

       do j = jsc, jec
          do i = isc, iec
             if( (wetv(i,j)>0.5).and.(miv(i,j)>MIV_MIN)) then ! timestep ice velocity (H&D eqn 22)
                rr       = cdw*dw*abs(cmplx(ui(i,j)-uo(i,j),vi(i,j)-vo(i,j)))*exp(sign(blturn*pi/180,cor(i,j))*(0.0,1.0))
                !
                ! first, timestep explicit parts (ice, wind & ocean part of water stress)
                !
                fxic_now = ( tmp1(i,j) + tmp2(i,j) + tmp3(i,j)*dxdy(i,j) - tmp4(i,j)*dydx(i,j) )/(dxv(i,j)*dyv(i,j)) 
                fyic_now = ( tmp6(i,j) + tmp7(i,j) + tmp3(i,j)*dydx (i,j) - tmp5(i,j)*dxdy(i,j))/(dxv(i,j)*dyv(i,j)) 

                ui(i,j) = ui(i,j)+(fxic_now+civ(i,j)*fxat(i,j)+ real(civ(i,j)*rr*cmplx(uo(i,j),vo(i,j))))*dtmiv(i,j)+sldx(i,j)
                vi(i,j) = vi(i,j)+(fyic_now+civ(i,j)*fyat(i,j)+aimag(civ(i,j)*rr*cmplx(uo(i,j),vo(i,j))))*dtmiv(i,j)+sldy(i,j)
                !
                ! second, timestep implicit parts (coriolis and ice part of water stress)
                !
                newuv = cmplx(ui(i,j),vi(i,j))/(1+dt_evp*(0.0,1.0)*cor(i,j)+civ(i,j)*rr*dtmiv(i,j))
                ui(i,j) = real(newuv); vi(i,j) = aimag(newuv)
                !
                ! sum for averages
                !
                fxic(i,j) = fxic(i,j) + fxic_now
                fyic(i,j) = fyic(i,j) + fyic_now
                fxoc(i,j) = fxoc(i,j) +  real(civ(i,j)*rr*cmplx(ui(i,j)-uo(i,j),vi(i,j)-vo(i,j)))
                fyoc(i,j) = fyoc(i,j) + aimag(civ(i,j)*rr*cmplx(ui(i,j)-uo(i,j),vi(i,j)-vo(i,j)))
                fxco(i,j) = fxco(i,j) - miv(i,j)*real ((0.0,1.0)*cor(i,j)*cmplx(ui(i,j),vi(i,j)))
                fyco(i,j) = fyco(i,j) - miv(i,j)*aimag((0.0,1.0)*cor(i,j)*cmplx(ui(i,j),vi(i,j)))              
             endif
          enddo
       enddo
    enddo
    !
    ! make averages
    !
    do j = jsc, jec
       do i = isc, iec
          if( (wetv(i,j)>0.5) .and. miv(i,j)>MIV_MIN ) then
             fxoc(i,j) = fxoc(i,j)/evp_sub_steps;  fyoc(i,j) = fyoc(i,j)/evp_sub_steps
             fxic(i,j) = fxic(i,j)/evp_sub_steps;  fyic(i,j) = fyic(i,j)/evp_sub_steps
             fxco(i,j) = fxco(i,j)/evp_sub_steps;  fyco(i,j) = fyco(i,j)/evp_sub_steps             
          endif
       enddo
    enddo

  end subroutine ice_dynamics

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! strain_angle                                                                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function strain_angle(ui, vi)
    real, dimension(isd:ied,jsd:jed), intent(in) :: ui, vi
    real, dimension(isc:iec,jsc:jec)             :: strn11, strn22, strn12, strain_angle

    integer :: i, j

    call set_strn(ui, vi, strn11, strn22, strn12)

    do j = jsc, jec
       do i = isc, iec
          if(strn11(i,j) + strn22(i,j) == 0.0 ) then
             strain_angle(i,j) = pi
          else
             strain_angle(i,j) = atan(((strn11(i,j)-strn22(i,j))**2+4*strn12(i,j)**2)**0.5/(strn11(i,j)+strn22(i,j)))
          endif
          if(strain_angle(i,j) < 0) then
             strain_angle(i,j) = 180 + strain_angle(i,j)*180/pi
          else
             strain_angle(i,j) = strain_angle(i,j)*180/pi
          endif
       enddo
    enddo

    return
  end function strain_angle

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! sigI - first stress invariant                                                !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function sigI(hi, ci, sig11, sig22, sig12)
    real, dimension(isc:iec,jsc:jec), intent(in) :: hi, ci, sig11, sig22, sig12
    real, dimension(isc:iec,jsc:jec)             :: sigI

    integer :: i, j

    sigI = ice_strength(hi,ci)

    do j = jsc, jec
       do i = isc, iec
          if(sigI(i,j) > 0.0) sigI(i,j) = (sig11(i,j) + sig22(i,j))/sigI(i,j)
       enddo
    enddo

    return
  end function sigI

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! sigII - second stress invariant                                              !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function sigII(hi, ci, sig11, sig22, sig12)
    real, dimension(isc:iec,jsc:jec), intent(in) :: hi, ci, sig11, sig22, sig12
    real, dimension(isc:iec,jsc:jec)             :: sigII

    integer :: i, j

    sigII = ice_strength(hi,ci)

    do j = jsc, jec
       do i = isc, iec
          if(sigII(i,j) > 0.0) sigII(i,j) = (((sig11(i,j)-sig22(i,j))**2+4*sig12(i,j)*sig12(i,j))/(sigII(i,j)**2))**0.5
       enddo
    enddo

    return
  end function sigII

end module ice_dyn_mod
