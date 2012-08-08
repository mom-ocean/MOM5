!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       THREE-LAYER VERTICAL THERMODYNAMICS                    !
!                                                                              !
! Reference:  M. Winton, 2000: "A reformulated three-layer sea ice model",     !
!            Journal of Atmospheric and Oceanic Technology, 17, 525-531.       !
!                                                                              !
!                                                                              !
!        -> +---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!       /   |         |                                                        !
!     hs    |  snow   | <- 0-heat capacity snow layer                          !
!       \   |         |                                                        !
!        => +---------+                                                        !
!       /   |         |                                                        !
!      /    |         | <- t1 - upper 1/2 ice temperature; this layer has      !
!     /     |         |         a variable (T/S dependent) heat capacity       !
!   hi      |...ice...|                                                        !
!     \     |         |                                                        !
!      \    |         | <- t2 - lower 1/2 ice temp. (fixed heat capacity)      !
!       \   |         |                                                        !
!        -> +---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                         Mike Winton (Michael.Winton)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! TK mod:  SLAB_ICE treatment modified to follow supersource
!          (after Bryan 1969).  The conductive heat
!           flux from ice to atmosphere is computed based on 
!           an effective ice thickness which ensures a minimum 
!           thickness of 1.7cm for the calculation.

module ice_thm_mod

use constants_mod, only : LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)

implicit none
private
public :: DS, DI, DW, MU_TS, TFI, ice_optics, ice3lay_temp, ice3lay_resize, &
          thm_pack, thm_unpack, ice_thm_param, e_to_melt
!
! properties of ice, snow, and seawater (NCAR CSM values)
!

real            :: KS    = 0.31      ! conductivity of snow - 0.31 W/(mK)
real, parameter :: DS    = 330.0     ! density of snow - 330 kg/(m^3)
real, parameter :: KI    = 2.03      ! conductivity of ice  - 2.03 W/(mK)
real, parameter :: DI    = 905.0     ! density of ice  - 905 kg/(m^3)
real, parameter :: CI    = 21e2      ! heat cap. of fresh ice - 2100 J/(kg K)
real, parameter :: SI    = 1.0       ! salinity of sea ice
real, parameter :: MU_TS = 0.054     ! relates freezing temp. to salinity
real, parameter :: TFI   = -MU_TS*SI ! sea ice freezing temp. = -mu*salinity
real, parameter :: CW    = 4.2e3     ! heat capacity of seawater
real, parameter :: DW    = 1030.0    ! density of water for waterline - kg/(m^3)

! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation
real            :: ALB_SNO = 0.85       ! albedo of snow (not melting)
real            :: ALB_ICE = 0.5826     ! albedo of ice (not melting)
real            :: PEN_ICE = 0.3        ! ice surface penetrating solar fraction
real            :: OPT_DEP_ICE = 0.67   ! ice optical depth (m)
real            :: T_RANGE_MELT = 1.0   ! melt albedos scaled in below melting T

real            :: H_LO_LIM = 0.0       ! hi/hs lower limit for temp. calc.
logical         :: SLAB_ICE = .false.   ! should we do old style GFDL slab ice?
!
! slab ice specific parameters
!
real, parameter :: CRIT_THICKNESS       = 1.00   
real, parameter :: T_RANGE              = 10.0
real, parameter :: MIN_ICE_ALB          = 0.55   ! coupled model uses 0.55
real, parameter :: MAX_ICE_ALB          = 0.80
real, parameter :: ALB_OCEAN            = 0.10

!
logical         :: CM2_BUGS = .false.   ! keep cm2 bugs for reproducibility

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_thm_param - set ice thermodynamic parameters                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_thm_param(alb_sno_in, alb_ice_in, pen_ice_in, opt_dep_ice_in, &
                         slab_ice_in, t_range_melt_in, cm2_bugs_in, ks_in, &
                         h_lo_lim_in                                          )
    real, intent(in)      :: alb_sno_in, alb_ice_in, pen_ice_in 
    real, intent(in)      :: opt_dep_ice_in, t_range_melt_in
logical, intent(in)   :: slab_ice_in
logical, intent(in)   :: cm2_bugs_in
    real, intent(in)   :: ks_in
    real, intent(in)   :: h_lo_lim_in

  ALB_SNO     = alb_sno_in
  ALB_ICE     = alb_ice_in
  PEN_ICE     = pen_ice_in
  OPT_DEP_ICE = opt_dep_ice_in
  SLAB_ICE    = slab_ice_in
  T_RANGE_MELT = t_range_melt_in

  CM2_BUGS    = cm2_bugs_in
  KS          = ks_in
  H_LO_LIM    = h_lo_lim_in

end subroutine ice_thm_param

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_optics - set albedo, penetrating solar, and ice/snow transmissivity      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_optics(alb, pen, trn, hs, hi, ts, tfw)
real, intent(  out) :: alb ! ice surface albedo (0-1)
real, intent(  out) :: pen ! fraction of down solar penetrating the ice
real, intent(  out) :: trn ! ratio of down solar at bottom to top of ice
real, intent(in   ) :: hs  ! snow thickness (m-snow)
real, intent(in   ) :: hi  ! ice thickness (m-ice)
real, intent(in   ) :: ts  ! surface temperature
real, intent(in   ) :: tfw ! seawater freezing temperature
real :: as, ai, cs
real :: thick_ice_alb, tcrit, fh


  if (SLAB_ICE) then
    tcrit = tfw - T_RANGE
       if (ts <= tcrit) then
          thick_ice_alb = MAX_ICE_ALB
       else if (ts >= tfw) then
          thick_ice_alb = MIN_ICE_ALB
       else
      thick_ice_alb = MAX_ICE_ALB + (MIN_ICE_ALB-MAX_ICE_ALB)*(ts-tcrit)/T_RANGE
       endif
  
       if (hi >= crit_thickness) then
          alb = THICK_ICE_ALB
       else
          alb = ALB_OCEAN + (thick_ice_alb-ALB_OCEAN)*sqrt(hi/CRIT_THICKNESS)
       endif

    pen = 0.0
    trn = 0.0

!! check for ice albdeos out of range (0 to 1)
!      if (alb.lt.0.0 .or.alb.gt.1.0) then
!         print *,'ice_optics: albedo out of range, alb_in=',alb_in, 'alb=',alb
!         print *,'ts=',ts,  'tfw=',tfw, 'tcrit=',tcrit
!         print *,'hi=',hi,  'thick_ice_alb=',thick_ice_alb
!         print *,'ALB_OCEAN=',ALB_OCEAN
!         print *,'MIN_ICE_ALB=',MIN_ICE_ALB, 'MAX_ICE_ALB=',MAX_ICE_ALB
!         print *,'T_RANGE,=',T_RANGE, 'CRIT_THICKNESS=',CRIT_THICKNESS
!         stop
!      end if
  
    return
  endif


!! 2007/04/11 Fix for thin ice negative ice albedos from Mike Winton
!!            Move ai calculation after if test

  as = ALB_SNO; ai = ALB_ICE
  cs = hs/(hs+0.02)                        ! thin snow partially covers ice

  fh = min(atan(5.0*hi)/atan(5.0*0.5),1.0) ! use this form from CSIM4 to
                                           ! reduce albedo for thin ice
  if (CM2_BUGS) then
    ai = fh*ai+(1-fh)*0.06                 ! reduce albedo for thin ice
    if (ts+T_RANGE_MELT > TFI) then        ! reduce albedo for melting as in
                                           ! CSIM4 assuming 0.53/0.47 vis/ir
       as = as-0.1235*min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
       ai = ai-0.075 *min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
    endif
  else
    if (ts+T_RANGE_MELT > TFI) then        ! reduce albedo for melting as in
                                           ! CSIM4 assuming 0.53/0.47 vis/ir
       as = as-0.1235*min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
       ai = ai-0.075 *min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
    endif
    ai = fh*ai+(1-fh)*0.06                 ! reduce albedo for thin ice
  end if

  alb = cs*as+(1-cs)*ai
  pen = (1-cs)*PEN_ICE
  trn = exp(-hi/OPT_DEP_ICE);

!! check for ice albdeos out of range (0 to 1)
! if (alb.lt.0.0 .or. alb.gt.1.0) then
!    print *,'ice_optics: albedo out of range, alb=',alb
!    print *,'cs=',cs,  'as=',as, 'ai=',ai
!    print *,'ts=',ts,  'fh=',fh, 'hs=',hs, 'hi=',hi, 'tfw=',tfw
!    print *,'ALB_SNO=',ALB_SNO,  'ALB_ICE=',ALB_ICE, 'T_RANGE_MELT,=',T_RANGE_MELT, 'TFI=',TFI
!    stop
! end if

  return

end subroutine ice_optics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice3lay_temp - ice & snow temp. change [Winton (2000) section 2.a]           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice3lay_temp(hs, hi, t1, t2, ts, A, B, I, tfw, fb, dt, tmelt, bmelt)
implicit none

real, intent(in   ) :: hs    ! snow thickness (m)
real, intent(in   ) :: hi    ! ice thickness (m)
real, intent(  out) :: ts    ! surface temperature (deg-C)
real, intent(inout) :: t1    ! upper ice temperature (deg-C)
real, intent(inout) :: t2    ! lower ice temperature (deg-C)
real, intent(in   ) :: A     ! net surface heat flux (+ up) at ts=0 (W/m^2)
real, intent(in   ) :: B     ! d(sfc heat flux)/d(ts) [W/(m^2 deg-C)]
real, intent(in   ) :: I     ! solar absorbed by upper ice (W/m^2)
real, intent(in   ) :: tfw   ! seawater freezing temperature (deg-C)
real, intent(in   ) :: fb    ! heat flux from ocean to ice bottom (W/m^2)
real, intent(in   ) :: dt    ! timestep (sec)
real, intent(inout) :: tmelt ! accumulated top melting energy  (J/m^2)
real, intent(inout) :: bmelt ! accumulated bottom melting energy (J/m^2)
!
! variables for temperature calculation [see Winton (1999) section II.A.]
! note:  here equations are multiplied by hi to improve thin ice accuracy
!
    real :: tsf, k12, hi2, a10, b10, a1, b1, c1
! TK Mods:
real :: hi_effective, hie
real :: KI_over_eps = 1.7065e-2     ! 5/2.93 from Bryan (1969); 
!                                 Value used in SS tsc.F (1.7065 cm) 
!                                  converted to meters...         

  if (SLAB_ICE) then
    hi_effective = hi + KI_over_eps     ! TK added
    ts = (KI*tfw-A*hi_effective)/(KI+B*hi_effective)     ! TK mod
    if (ts > 0.0) then       ! surface melting conditions
       ts = 0.0
       if (hi>0.0) tmelt = tmelt + (KI*tfw/hi_effective-A)*dt     ! TK mod
    endif
    if (hi>0.0) then
       bmelt = bmelt + (fb-KI*(tfw-ts)/hi_effective)*dt     ! TK mod
    else
       bmelt = bmelt + (fb-A-B*tfw)*dt
    endif
    return
  endif

  if (hs > 0.0) then
    TSF = 0.0
  else
    TSF = TFI
  endif

  hie = max(hi, H_LO_LIM); ! prevent thin ice inaccuracy (mw)

  !
  ! Compute upper ice and surface temperatures
  !
  K12 = 4*KI*KS/(KS+4*KI*hs/hie)
  hi2 = hie*hie

  A10 = DI*hi2*CI/(2*dt) + 2*KI*(4*dt*2*KI+DI*hi2*CI)/(6*dt*2*KI+DI*hi2*CI)
  B10 = -DI*hi2*(CI*t1+LI*TFI/t1)/(2*dt) - I*hie                       &
        -2*KI*(4*dt*2*KI*tfw+DI*hi2*CI*t2)/(6*dt*2*KI+DI*hi2*CI)

  A1 = A10+K12*B*hie/(K12+B*hie)
  B1 = B10+A*K12*hie/(K12+B*hie)
  C1  = DI*hi2*LI*TFI/(2*dt)
  t1 = -(sqrt(B1*B1-4*A1*C1)+B1)/(2*A1)
  ts = (K12*t1-A*hie)/(K12+B*hie)
   
  if (ts > tsf) then       ! slightly different equation for melting conditions
    A1 = A10+K12
    B1 = B10-K12*tsf
    t1 = -(sqrt(B1*B1-4*A1*C1)+B1)/(2*A1)
    ts = tsf
    tmelt = tmelt + (K12*(t1-ts)/hie-(A+B*ts))*dt
  endif
  !
  ! set lower ice temp. -- use tfw as reference for thin ice precision
  !
  t1 = t1-tfw; t2 = t2-tfw;
  t2 = (2*dt*2*KI*t1+DI*hi2*CI*t2)/(6*dt*2*KI+DI*hi2*CI)
  t1 = t1+tfw; t2 = t2+tfw;

  bmelt = bmelt + (fb+4*KI*(t2-tfw)/hie)*dt

  if (t2 > TFI) then ! put excess lower ice energy into bmelt
    bmelt = bmelt + e_to_melt(h2=hie/2,t2=TFI) - e_to_melt(h2=hie/2,t2=t2)
    t2 = TFI
  endif

  if (t1 > TFI) then ! put excess upper ice energy into tmelt
    tmelt = tmelt + e_to_melt(h1=hie/2,t1=TFI) - e_to_melt(h1=hie/2,t1=t1)
    t1 = TFI
  endif

! if (tmelt<0) then
!    print *,'neg. tmelt=',tmelt,ts,t1,hs,hi,K12*(t1-ts)/hi,-(A+B*ts)
!    print *,'K12=',K12
!    print *,'A/B/I=',A,B,I
!    print *,'A/B/C=',A1,B1,C1,A1*t1*t1+B1*t1+C1,A1*t1*t1,B1*t1,C1
!    stop
! end if
! call thm_checkout(ts, hs, hi, t1, t2, bmelt, tmelt)

  return
end subroutine ice3lay_temp

subroutine thm_checkout(ts, hs, hi, t1, t2, bmelt, tmelt)
    real, intent(in) :: ts, hs, hi, t1, t2, bmelt, tmelt

  if (ts>0 .or. t1>TFI .or. t2>0 .or. hs<0 .or. hs>100 .or. hi<0 .or. hi>100 &
           .or. abs(bmelt)>100*DI*LI .or. tmelt<0 .or. tmelt>100*DI*LI) then
      print *,'UNREASONABLE ICE: hs=',hs,'hi=',hi,'t1=',t1,'t2=',t2,'ts=', &
              ts,'tmelt=',tmelt,'bmelt=',bmelt
  end if
end subroutine thm_checkout

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! e_to_melt - energy needed to melt a given snow/ice configuration             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function e_to_melt(hs, h1, t1, h2, t2)
    real, intent(in), optional :: hs, h1, t1, h2, t2
    real                       :: e_to_melt

  e_to_melt = 0.0
  if (present(hs))              e_to_melt = e_to_melt+DS*LI*hs
  if (present(h1).and.present(t1)) then
                                e_to_melt = e_to_melt+DI*h1*(CI-LI/t1)*(TFI-t1)
  endif
  if (present(h2).and.present(t2)) then
                                e_to_melt = e_to_melt+DI*h2*(LI+CI*(TFI-t2))
  endif
end function e_to_melt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! add_to_top - add some ice to the top ice layer                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine add_to_top(h, t, h1, t1)
    real,    intent(in) :: h, t
    real, intent(inout) :: h1, t1
    real                :: f1

  f1 = h1/(h1+h)
  t1 = f1*(t1+LI*TFI/(CI*t1))+(1-f1)*t
  t1 = (t1-sqrt(t1*t1-4*TFI*LI/CI))/2
  h1 = h1+h
  return
end subroutine add_to_top

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! add_to_bot - add some ice to the bottom ice layer                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine add_to_bot(h, t, h2, t2)
    real,    intent(in) :: h, t
    real, intent(inout) ::  h2, t2

  t2 = (h2*t2+h*t)/(h2+h)
  h2 = h2+h
  return
end subroutine add_to_bot

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! even_up - transfer mass/energy between ice layers to maintain equal thickness!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine even_up(h1, t1, h2, t2)
    real,intent(inout) :: h1, t1, h2, t2
    real               :: dh

  if (h1 > (h1+h2)/2) then
    call add_to_bot(h1-(h1+h2)/2, t1+LI*TFI/(CI*t1), h2, t2)
    h1 = h2
  else if (h2 > (h1+h2)/2) then
    call add_to_top(h2-(h1+h2)/2, t2, h1, t1)
    h2 = h1
  endif
  if (t2>TFI) then                 ! use extra energy to melt both layers evenly
    dh = h2*CI*(t2-TFI)*t1/(LI*t1+(CI*t1-LI)*(TFI-t1))
    t2 = TFI
    h1 = h1-dh
    h2 = h2-dh
  endif
  return
end subroutine even_up

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice3lay_resize - ice & snow thickness change [Winton (1998) section II.B.]   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice3lay_resize(hs, hi, t1, t2, snow, frazil, evap, tmelt, bmelt, &
                          tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn,       &
                          snow_to_ice, bablt                                )
real, intent(inout) :: hs          ! snow thickness (m-snow)
real, intent(inout) :: hi          ! ice thickness (m-ice)
real, intent(inout) :: t1          ! temperature of upper ice (deg-C)
real, intent(inout) :: t2          ! temperature of lower ice (deg-C)
real, intent(in   ) :: snow        ! new snow (kg/m^2-snow)
real, intent(in   ) :: frazil      ! frazil in energy units
real, intent(in   ) :: evap        ! ice evaporation (kg/m^2)
real, intent(in   ) :: tmelt       ! top melting energy (J/m^2)
real, intent(in   ) :: bmelt       ! bottom melting energy (J/m^2)
real, intent(in   ) :: tfw         ! seawater freezing temperature (deg-C)
real, intent(  out) :: heat_to_ocn ! energy left after ice all melted (J/m^2)
real, intent(  out) :: h2o_to_ocn  ! liquid water flux to ocean (kg/m^2)
real, intent(  out) :: h2o_from_ocn! evaporation flux from ocean (kg/m^2)
real, intent(  out) :: snow_to_ice ! snow below waterline becomes ice
real, intent(  out), optional :: bablt ! bottom ablation (kg/m^2)

real h1, h2, dh, hw

  heat_to_ocn  = 0.0
  h2o_to_ocn   = DS*hs+DI*hi+snow-evap ! - from ice at end gives ocean h2o flux
  h2o_from_ocn = 0.0
  snow_to_ice  = 0.0

  if (SLAB_ICE) then
    !
    ! add snow and frazil
    !
    hi = hi + snow/DI + frazil/(DI*LI)
    t1 = tfw; t2 = tfw;
    !
    ! atmospheric evaporation
    !
    if (evap <= hi*DI) then
      hi = hi - evap/DI
    else
      h2o_from_ocn = evap-hi*DI
      hi = 0.0
    end if
    !
    ! ... melting
    !
    hi = hi - (tmelt+bmelt)/(DI*LI)
    if (hi<0.0) then
       heat_to_ocn = -hi*DI*LI
       hi = 0.0
    else
       heat_to_ocn = 0.0
    endif
  
    h2o_to_ocn = h2o_to_ocn+h2o_from_ocn ! reset mark for leftover evap thru ice
    h2o_to_ocn = h2o_to_ocn-DI*hi-DS*hs  ! hs should be zero
    if (present(bablt)) bablt = bmelt/LI
    return
  endif

  h1 = hi/2
  h2 = hi/2
  !
  ! add snow ...
  !
  hs = hs + snow/DS
  !
  ! ... and frazil
  !
  call add_to_bot(frazil/e_to_melt(h2=1.0,t2=tfw), tfw, h2, t2)
  !
  ! atmospheric evaporation
  !
  if (evap <= hs*DS) then
    hs = hs - evap/DS
  else if (evap-hs*DS<=h1*DI) then
    hs = 0.0
    h1 = h1-(evap-DS*hs)/DI
  else if (evap-hs*DS-h1*DI<=h2*DI) then
    hs = 0.0
    h1 = 0.0
    h2 = h2 - (evap-hs*DS-h1*DI)/DI
  else
    h2o_from_ocn = evap-hs*DS-(h1+h2)*DI
    hs = 0.0
    h1 = 0.0
    h2 = 0.0
  end if

  if (bmelt < 0.0) call add_to_bot(-bmelt/e_to_melt(h2=1.0,t2=tfw), tfw, h2, t2)

  if (h1 == 0.0) t1 = tfw  ! need this, below we divide by t1 even when h1 == 0

  !
  ! apply energy fluxes ... top ...
  !
  if (tmelt <= e_to_melt(hs)) then
    hs = hs - tmelt/e_to_melt(hs=1.0)
  else if (tmelt <= e_to_melt(hs,h1,t1)) then
    h1 = h1 - (tmelt-e_to_melt(hs))/e_to_melt(h1=1.0,t1=t1)
    hs = 0.0
  else if (tmelt <= e_to_melt(hs,h1,t1,h2,t2)) then
    h2 = h2 - (tmelt-e_to_melt(hs,h1,t1))/e_to_melt(h2=1.0,t2=t2)
    hs = 0.0
    h1 = 0.0
  else
   heat_to_ocn = heat_to_ocn+tmelt-e_to_melt(hs,h1,t1,h2,t2)
   hs = 0.0
   h1 = 0.0
   h2 = 0.0
  endif
  !
  ! ... and bottom
  !
  if (present(bablt)) bablt = DS*hs+DI*(h1+h2)
  if (bmelt > 0.0) then
    if (bmelt < e_to_melt(h2=h2,t2=t2)) then
      h2 = h2 - bmelt/e_to_melt(h2=1.0,t2=t2)
    else if (bmelt < e_to_melt(h1=h1,t1=t1,h2=h2,t2=t2)) then
      h1 = h1-(bmelt-e_to_melt(h2=h2,t2=t2))/e_to_melt(h1=1.0,t1=t1)
      h2 = 0.0
    else if (bmelt < e_to_melt(hs,h1,t1,h2,t2)) then
      hs = hs - (bmelt-e_to_melt(h1=h1,t1=t1,h2=h2,t2=t2))/e_to_melt(hs=1.0)
      h1 = 0.0
      h2 = 0.0
    else
      heat_to_ocn = heat_to_ocn+bmelt-e_to_melt(hs,h1,t1,h2,t2)
      hs = 0.0
      h1 = 0.0
      h2 = 0.0
    endif
  endif
  if (present(bablt)) bablt = bablt-DS*hs-DI*(h1+h2)

  hi = h1 + h2
  hw = (DI*hi+DS*hs)/DW
  if (hw>hi) then           ! convert snow to ice to maintain ice at waterline
    snow_to_ice = (hw-hi)*DI
    hs = hs - snow_to_ice/DS
    call add_to_top(hw-hi, TFI, h1, t1)
  endif

  call even_up(h1, t1, h2, t2)
  hi = h1+h2
  if (hi==0.0) then
    t1 = 0.0
    t2 = 0.0
  endif

  h2o_to_ocn = h2o_to_ocn+h2o_from_ocn ! correct mark for leftover evap thru ice
  h2o_to_ocn = h2o_to_ocn-DS*hs-DI*hi

  return
end subroutine ice3lay_resize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! thm_pack - form conserved quantities for two-dimensional advection           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine thm_pack(cn, hs, hi, t1, t2)
real, dimension(:,:,:), intent(inout) :: cn, hs, hi, t1, t2

    integer :: i, j, k
    real    :: tmp

    tmp = MU_TS*SI*LI/CI

    do k = 1, size(hi,3)
       do j = 1, size(hi,2)
          do i = 1, size(hi,1)
             if (hi(i,j,k)>0.0) then
                hi(i,j,k) = cn(i,j,k)*hi(i,j,k)
                hs(i,j,k) = cn(i,j,k)*hs(i,j,k)
                t1(i,j,k) = (t1(i,j,k)-tmp/t1(i,j,k))*hi(i,j,k)
                t2(i,j,k) = t2(i,j,k)*hi(i,j,k)
             else 
                cn(i,j,k) = 0.0
                hi(i,j,k) = 0.0
                hs(i,j,k) = 0.0
                t1(i,j,k) = 0.0
                t2(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

  return

end subroutine thm_pack

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! thm_unpack - reform ice properties from advectively conserved quantities     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine thm_unpack(cn, hs, hi, t1, t2)
real, dimension(:,:,:), intent(inout) :: cn, hs, hi, t1, t2
    integer :: i, j, k
    real    :: tmp

    tmp = 4*MU_TS*SI*LI/CI

    do k = 1, size(hi,3)
       do j = 1, size(hi,2)
          do i = 1, size(hi,1)
             if (hi(i,j,k)>0.0) then
                t1(i,j,k) = t1(i,j,k)/hi(i,j,k)
                t1(i,j,k) = 0.5*(t1(i,j,k)-sqrt(t1(i,j,k)*t1(i,j,k)+tmp))
                t2(i,j,k) = t2(i,j,k)/hi(i,j,k)
                hi(i,j,k) = hi(i,j,k)/cn(i,j,k)
                hs(i,j,k) = hs(i,j,k)/cn(i,j,k)
             else
                cn(i,j,k) = 0.0
                hi(i,j,k) = 0.0
                hs(i,j,k) = 0.0
                t1(i,j,k) = 0.0
                t2(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

  return
end subroutine thm_unpack

end module ice_thm_mod
