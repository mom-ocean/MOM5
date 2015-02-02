      module MO_SETSOX_MOD

implicit none
character(len=128), parameter :: version     = '$Id: mo_setsox.F90,v 19.0 2012/01/06 20:34:08 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine setsox( press, plonl, dtime, tfld, qfld, &
                         lwc, xhnm, &
                         qin, retain_cm3_bugs )

!-----------------------------------------------------------------------      
!          ... Compute heterogeneous reactions of SOX
!
!       (0) using initial PH to calculate PH
!           (a) HENRY's law constants
!           (b) PARTIONING
!           (c) PH values
!
!       (1) using new PH to repeat
!           (a) HENRY's law constants
!           (b) PARTIONING
!           (c) REACTION rates
!           (d) PREDICTION
!-----------------------------------------------------------------------      

      use mo_chem_utls_mod, only : get_spc_ndx

      implicit none
!-----------------------------------------------------------------------      
!      ... Dummy arguments
!-----------------------------------------------------------------------      
      integer, intent(in)  ::    plonl               ! number of local longitude points
      real, intent(in)     ::    dtime               ! time step (sec)
      real, intent(inout)  ::    qin(:,:,:)          ! xported species ( vmr )
      real, intent(in)     ::    xhnm(:,:)           ! total atms density ( /cm**3)
      real, dimension(:,:), intent(in) ::  &
                               tfld, &               ! temperature
                               qfld, &               ! specific humidity( kg/kg )
                               lwc,  &               ! cloud liquid water content (kg/kg)
                               press                 ! midpoint pressure ( Pa )
      logical, intent(in) ::   retain_cm3_bugs       ! retain cm3 bugs ?

!-----------------------------------------------------------------------      
!      ... Local variables
!
!           xhno3 ... in mixing ratio
!-----------------------------------------------------------------------      
      integer, parameter :: itermax = 20
      real, parameter ::  ph0 = 5.0  ! INITIAL PH VALUES
      real, parameter ::  const0 = 1.e3/6.023e23
      real, parameter ::  xa0 = 11.,   &
                          xb0 = -.1,   &
                          xa1 = 1.053, &
                          xb1 = -4.368,&
                          xa2 = 1.016, &
                          xb2 = -2.54, &
                          xa3 = .816e-32, &
                          xb3 = .259
      real, parameter ::  kh0 = 9.e3, &           ! HO2(g)          -> Ho2(a)
                          kh1 = 2.05e-5, &        ! HO2(a)          -> H+ + O2-
                          kh2 = 8.6e5,   &        ! HO2(a) + ho2(a) -> h2o2(a) + o2
                          kh3 = 1.e8,    &        ! HO2(a) + o2-    -> h2o2(a) + o2
                          Ra = 8314./101325., &   ! universal constant   (atm)/(M-K)
                          xkw = 1.e-14            ! water acidity
      real, parameter :: small_value = 1.e-20

      integer    ::      k, i, iter
      integer    ::      plev
      real       ::      wrk, delta
      real       ::      xph0, xk, xe, x2
      real       ::      tz, xl, px, qz, pz, es, qs, patm
      real       ::      Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
      real       ::      hno3g, nh3g, so2g, h2o2g, co2g, o3g
      real       ::      rah2o2, rao3, pso4, ccc
      real       ::      xx0, yy1, xkp
      real       ::      cnh3, chno3, com, com1, com2, xra
      real       ::      RH
      integer    ::      ox_ndx, hno3_ndx, h2o2_ndx, so2_ndx, so4_ndx, &
                         nh3_ndx, nh4no3_ndx, ho2_ndx

!-----------------------------------------------------------------------      
!            for Ho2(g) -> H2o2(a) formation 
!            schwartz JGR, 1984, 11589
!-----------------------------------------------------------------------      
      real       ::      kh4    ! kh2+kh3
      real       ::      xam    ! air density /cm3
      real       ::      ho2s   ! ho2s = ho2(a)+o2-
      real       ::      r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
      real       ::      r2h2o2 ! prod(h2o2) by ho2 in mix/s

      real, dimension(SIZE(tfld,1),SIZE(tfld,2))  :: &
                         xhno3, xh2o2, xso2, xso4,&
                         xnh3, xo3,         &
                         xlwc, cfact,       &
                         xph, xant, xho2,         &
                         hehno3, &            ! henry law const for hno3
                         heh2o2, &            ! henry law const for h2o2
                         heso2,  &            ! henry law const for so2
                         henh3,  &            ! henry law const for nh3
                         heo3                 ! henry law const for nh3
      real, dimension(plonl)  :: &
                         t_fac
      logical :: converged


!      ox_ndx = get_spc_ndx( 'OX' )
! for ox budget (jmao, 1/7/2011)
    if (retain_cm3_bugs) then
      ox_ndx = get_spc_ndx( 'OX' )
    else
      ox_ndx = get_spc_ndx( 'O3' )
    endif
      hno3_ndx = get_spc_ndx( 'HNO3' )
      h2o2_ndx = get_spc_ndx( 'H2O2' )
      so2_ndx = get_spc_ndx( 'SO2' )
      so4_ndx = get_spc_ndx( 'SO4' )
      nh3_ndx = get_spc_ndx( 'NH3' )
      nh4no3_ndx = get_spc_ndx( 'NH4NO3' )
      ho2_ndx = get_spc_ndx( 'HO2' )

      plev = SIZE(tfld,2)
      
!-----------------------------------------------------------------
!       ... NOTE: The press array is in pascals and must be
!                 mutiplied by 10 to yield dynes/cm**2.
!-----------------------------------------------------------------

!==================================================================
!       ... First set the PH
!==================================================================
!      ... Initial values
!           The values of so2, so4 are after (1) SLT, and CHEM
!-----------------------------------------------------------------
      xph0 = 10.**(-ph0)                      ! initial PH value
      do k = 1,plev
!        precip(:,k) = cmfdqr(:,k) + rain(:,k) - evapr(:,k)
      end do

      do k = 1,plev
          cfact(1:,k) = xhnm(1:,k) &             ! /cm3(a)  
                            * 1.e6          &             ! /m3(a)
                            * 1.38e-23/287. &             ! Kg(a)/m3(a)
                            * 1.e-3                       ! Kg(a)/L(a)
      end do

      do k = 1,plev
         xph(:,k) = xph0                                    ! initial PH value
         xlwc(:,k) = lwc(:,k) *cfact(:,k)           ! cloud water  L(water)/L(air)
!        xrain(:,k) = rain(:,k) *cfact(:,k)         ! RAIN  water  L(water)/L(air)
         if( hno3_ndx > 0 ) then
            xhno3(:,k) = qin(:,k,hno3_ndx)                  ! mixing ratio
         else
            xhno3(:,k) = 0.
         end if
         if( h2o2_ndx > 0 ) then
            xh2o2(:,k) = qin(:,k,h2o2_ndx)                  ! mixing ratio
         else
            xh2o2(:,k) = 0.
         end if
         if( so2_ndx > 0 ) then
            xso2(:,k) = qin(:,k,so2_ndx)                   ! mixing ratio
         else
            xso2(:,k) = 0.
         end if
         if( so4_ndx > 0 ) then
            xso4(:,k) = qin(:,k,so4_ndx)                   ! mixing ratio
         else
            xso4(:,k) = 0.
         end if
         if( nh3_ndx > 0 ) then
            xnh3(:,k) = qin(:,k,nh3_ndx)                   ! mixing ratio
         else
            xnh3(:,k) = 0.
         end if
         if( nh4no3_ndx > 0 ) then
            xant(:,k) = qin(:,k,nh4no3_ndx)                   ! mixing ratio
         else
            xant(:,k) = 0.
         end if
         if( ox_ndx > 0 ) then
            xo3(:,k) = qin(:,k,ox_ndx)                    ! mixing ratio
         else
            xo3(:,k) = 0.
         end if
         if( ho2_ndx > 0 ) then
            xho2(:,k) = qin(:,k,ho2_ndx)                   ! mixing ratio
         else
            xho2(:,k) = 0.
         end if
      end do 

!-----------------------------------------------------------------
!       ... Temperature dependent Henry constants
!-----------------------------------------------------------------
      do k = 1,plev                                             !! plev loop for STEP 0
         do i = 1,plonl
            xl = xlwc(i,k) 
            if( xl >= 1.e-8 ) then
               t_fac(i) = 1. / tfld(i,k) - 1. / 298.
!-----------------------------------------------------------------------      
!        ... hno3
!-----------------------------------------------------------------------      
               do iter = 1,itermax
                  xk = 2.1e5 *EXP( 8700.*t_fac(i) )
                  xe = 15.4
                  hehno3(i,k)  = xk*(1. + xe/xph(i,k))
!-----------------------------------------------------------------------      
!         ... h2o2
!-----------------------------------------------------------------------      
                  xk = 7.4e4   *EXP( 6621.*t_fac(i) )
                  xe = 2.2e-12 *EXP(-3730.*t_fac(i) )
                  heh2o2(i,k)  = xk*(1. + xe/xph(i,k))
!-----------------------------------------------------------------------      
!          ... so2
!-----------------------------------------------------------------------      
                  xk = 1.23  *EXP( 3120.*t_fac(i) )
                  xe = 1.7e-2*EXP( 2090.*t_fac(i) )
                  x2 = 6.0e-8*EXP( 1120.*t_fac(i) )
                  wrk = xe/xph(i,k)
                  heso2(i,k)  = xk*(1. + wrk*(1. + x2/xph(i,k)))
!-----------------------------------------------------------------------      
!          ... nh3
!-----------------------------------------------------------------------      
                  xk = 58.   *EXP( 4085.*t_fac(i) )
                  xe = 1.7e-5*EXP(-4325.*t_fac(i) )
                  henh3(i,k)  = xk*(1. + xe*xph(i,k)/xkw)
!-----------------------------------------------------------------
!       ... Partioning and effect of pH 
!-----------------------------------------------------------------
                  pz = .01*press(i,k)       !! pressure in mb
                  tz = tfld(i,k)
                  patm = pz/1013.
                  xam  = press(i,k)/(1.38e-23*tz)  !air density /M3
!-----------------------------------------------------------------
!        ... hno3
!-----------------------------------------------------------------
                  px = hehno3(i,k) * Ra * tz * xl
                  hno3g = xhno3(i,k)/(1. + px)
                  xk = 2.1e5 *EXP( 8700.*t_fac(i) )
                  xe = 15.4
                  Ehno3 = xk*xe*hno3g *patm
!-----------------------------------------------------------------
!          ... so2
!-----------------------------------------------------------------
                  px = heso2(i,k) * Ra * tz * xl
                  so2g =  xso2(i,k)/(1.+ px)
                  xk = 1.23  *EXP( 3120.*t_fac(i) )
                  xe = 1.7e-2*EXP( 2090.*t_fac(i) )
                  Eso2 = xk*xe*so2g *patm
!-----------------------------------------------------------------
!          ... nh3
!-----------------------------------------------------------------
                  px = henh3(i,k) * Ra * tz * xl
                  nh3g = xnh3(i,k)/(1.+ px)
                  xk = 58.   *EXP( 4085.*t_fac(i) )
                  xe = 1.7e-5*EXP( -4325.*t_fac(i) )
                  Enh3 = xk*xe*nh3g/xkw *patm
!-----------------------------------------------------------------
!        ... h2o effects
!-----------------------------------------------------------------
                  Eh2o = xkw
!-----------------------------------------------------------------
!        ... co2 effects
!-----------------------------------------------------------------
                  co2g = 330.e-6                            !330 ppm = 330.e-6 atm
                  xk = 3.1e-2*EXP( 2423.*t_fac(i) )
                  xe = 4.3e-7*EXP(-913. *t_fac(i) )
                  Eco2 = xk*xe*co2g  *patm
!-----------------------------------------------------------------
!        ... PH cal
!-----------------------------------------------------------------
                  com2 = (Eh2o + Ehno3 + Eso2 + Eco2)  &
                       / (1. + Enh3 )
                  com2 = MAX( com2,1.e-20 )
                  xph(i,k) = SQRT( com2 )
!-----------------------------------------------------------------
!         ... Add so4 effect
!-----------------------------------------------------------------
                  Eso4 = xso4(i,k)*xhnm(i,k)   &         ! /cm3(a)
                        *const0/xl
                  xph(i,k) =  MIN( 1.e-2,MAX( 1.e-7,xph(i,k) + 2.*Eso4 ) )
                  if( iter > 1 ) then
                     if ( ABS(delta) > 1.e-40 ) then
                        delta = ABS( (xph(i,k) - delta)/delta )
                     else
                        delta = 0.
                     end if
                     converged = delta < .01
                     if( converged ) then
                        exit
                     else
                        delta = xph(i,k)
                     end if
                  else
                     delta = xph(i,k)
                  end if
               end do
               if( .not. converged ) then
                  write(*,*) 'SETSOX: pH failed to converge @ (',i,',',k,'), % change=', &
                              100.*delta
               end if
            else
               xph(i,k) =  1.e-7
            end if
         end do
      end do  ! end plev loop for STEP 0
!     do file = 1,match_file_cnt
!     call OUTFLD( 'PH', xph,  plonl, ip, lat, file )
!     ENDDO

!==============================================================
!          ... Now use the actual PH
!==============================================================
      do k = 1,plev
         do i = 1,plonl
            t_fac(i) = 1. / tfld(i,k) - 1. / 298.
            tz = tfld(i,k)
            xl = xlwc(i,k)
            patm = press(i,k)/101300.        ! press is in pascal
            xam  = press(i,k)/(1.38e-23*tz)  ! air density /M3

!-----------------------------------------------------------------
!         ... hno3
!-----------------------------------------------------------------
            xk = 2.1e5 *EXP( 8700.*t_fac(i) )
            xe = 15.4
            hehno3(i,k)  = xk*(1. + xe/xph(i,k))

!-----------------------------------------------------------------
!        ... h2o2
!-----------------------------------------------------------------
            xk = 7.4e4   *EXP( 6621.*t_fac(i) )
            xe = 2.2e-12 *EXP(-3730.*t_fac(i) )
            heh2o2(i,k)  = xk*(1. + xe/xph(i,k))

!-----------------------------------------------------------------
!         ... so2
!-----------------------------------------------------------------
            xk = 1.23  *EXP( 3120.*t_fac(i) )
            xe = 1.7e-2*EXP( 2090.*t_fac(i) )
            x2 = 6.0e-8*EXP( 1120.*t_fac(i) )

            wrk = xe/xph(i,k)
            heso2(i,k)  = 1.e2  !xk*(1. + wrk*(1. + x2/xph(i,k)))

!-----------------------------------------------------------------
!          ... nh3
!-----------------------------------------------------------------
            xk = 58.   *EXP( 4085.*t_fac(i) )
            xe = 1.7e-5*EXP(-4325.*t_fac(i) )
            henh3(i,k)  = xk*(1. + xe*xph(i,k)/xkw)

!-----------------------------------------------------------------
!        ... o3
!-----------------------------------------------------------------
            xk = 1.15e-2 *EXP( 2560.*t_fac(i) )
            heo3(i,k) = xk

!------------------------------------------------------------------------
!       ... for Ho2(g) -> H2o2(a) formation 
!           schwartz JGR, 1984, 11589
!------------------------------------------------------------------------
            kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1. + kh1/xph(i,k))**2)
            ho2s = kh0*xho2(i,k)*patm*(1. + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
            r1h2o2 = kh4*ho2s*ho2s                         ! prod(h2o2) in mole/L(w)/s
            r2h2o2 = r1h2o2*xlwc(i,k)  &                   ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                           *const0     &                   ! mole/fm3(a)/s * 1.e-3       = mole/cm3(a)/s
                           /xam                            ! /cm3(a)/s    / air-den     = mix-ratio/s
            xh2o2(i,k) = xh2o2(i,k) + r2h2o2*dtime         ! updated h2o2 by het production

!-----------------------------------------------
!       ... Partioning 
!-----------------------------------------------
!------------------------------------------------------------------------
!        ... h2o2
!------------------------------------------------------------------------
            px = heh2o2(i,k) * Ra * tz * xl
            h2o2g =  xh2o2(i,k)/(1.+ px)

!------------------------------------------------------------------------
!         ... so2
!------------------------------------------------------------------------
            px = heso2(i,k) * Ra * tz * xl
            so2g =  xso2(i,k)/(1.+ px)

!------------------------------------------------------------------------
!         ... o3 ============
!------------------------------------------------------------------------
            px = heo3(i,k) * Ra * tz * xl
            o3g =  xo3(i,k)/(1.+ px)

!-----------------------------------------------
!       ... Aqueous phase reaction rates
!           SO2 + H2O2 -> SO4
!           SO2 + O3   -> SO4
!-----------------------------------------------
          
!------------------------------------------------------------------------
!       ... S(IV) (HSO3) + H2O2
!------------------------------------------------------------------------
            rah2o2 = 8.e4 * EXP( -3650.*t_fac(i) )  &
                   / (.1 + xph(i,k))

!------------------------------------------------------------------------
!        ... S(IV)+ O3
!------------------------------------------------------------------------
            rao3   = 4.39e11 * EXP(-4131./tz)  &
                  + 2.56e3  * EXP(-996. /tz) /xph(i,k)

!-----------------------------------------------------------------
!       ... Prediction after aqueous phase
!       so4
!       When Cloud is present 
!   
!       S(IV) + H2O2 = S(VI)
!       S(IV) + O3   = S(VI)
!
!       reference:
!           (1) Seinfeld
!           (2) Benkovitz
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... S(IV) + H2O2 = S(VI)
!-----------------------------------------------------------------
            if( xl >= 1.e-8 ) then                          ! when cloud is present
               pso4 = rah2o2 * heh2o2(i,k)*h2o2g  &
                             * heso2(i,k) *so2g             ! [M/s]
               pso4 = pso4       &                          ! [M/s] =  [mole/L(w)/s]
                    * xlwc(i,k)  &                          ! [mole/L(a)/s]
                    / const0     &                          ! [/L(a)/s]
                    / xhnm(i,k)                             ! [mixing ratio/s]

          ccc = pso4*dtime


               ccc = MAX( MIN( ccc, xso2(i,k), xh2o2(i,k) ), 0. )
                  xso4(i,k)  = xso4(i,k)  + ccc
               xh2o2(i,k) = MAX( xh2o2(i,k) - ccc, small_value )
               xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )

!-----------------------------------------------
!       ... S(IV) + O3 = S(VI)
!-----------------------------------------------
               pso4 = rao3 * heo3(i,k)*o3g * heso2(i,k)*so2g       ! [M/s]
               pso4 = pso4        &                                ! [M/s] =  [mole/L(w)/s]
                    * xlwc(i,k)   &                                ! [mole/L(a)/s]
                    / const0      &                                ! [/L(a)/s]
                    / xhnm(i,k)                                    ! [mixing ratio/s]

          ccc = pso4*dtime

               ccc = MAX( MIN( ccc, xso2(i,k) ), 0. )
                  xso4(i,k)  = xso4(i,k)  + ccc
               xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )
            end if                                               ! when cloud is present

!-----------------------------------------------------------------
!       ... Formation of NH4+ + SO4=
!           to balance 1 SO4= should take 2 NH4+
!           According to Dentener and Crutzen (1994) JAC 331
!           the neutralization of sulfuric acid by NH3
!           is (NH4)1.5 H0.5(SO4)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... Formation of AMMONIUM NITRATE (ANT)
!           Calculate reaction coefficient NH3(g)+HNO3(g)=NH4(+)NO3(-) 
!                                                   Kp(ppb**2)
!      * Kp is calculated according to
!        Stelson and Seinfeld Atm. Env 16 983, 1982
!        Seinfeld (1986) 
!-----------------------------------------------------------------
            qz = qfld(i,k)             ! H2O mass mxing ratio Kg/Kg
            pz = .01*press(i,k)        ! pressure in mb
 
!-----------------------------------------------------------------
!        ... Calculate RH
!-----------------------------------------------------------------
            wrk = tz - 273.
            es = 6.11*10.**(7.63*wrk/(241.9 + wrk))            ! Magnus EQ
            qs = .622*es/pz                                    ! sat mass mix (H2O)
            RH = 100.*qz/qs                                    ! relative huminity(%)
            RH = MIN( 100.,MAX( RH,0. ) )
  
            xx0 = xa0 + xb0*RH

            if( RH >= 90. ) then
               yy1 = xa1*EXP( xb1/xx0 )
            else
               yy1 = xa2*EXP( xb2/xx0 )
            end if            

            xkp = yy1*(xa3*EXP( xb3*tz )/.7) &    ! ppb**2
                    * 1.e-18                      ! mixing ratio

            cnh3 = xnh3(i,k)
            chno3 = xhno3(i,k)
            com = cnh3*chno3

            com1 = (cnh3 + chno3)**2 - 4.*(cnh3*chno3 - xkp)
            com1 = MAX( com1,1.e-30 )

            if( com >= xkp ) then   ! NH4NO3 is formed
               xra = .5*(cnh3 + chno3 - SQRT(com1))
!-----------------------------------------------------------------
!        ... xra =0.0 for not forming ANT
!-----------------------------------------------------------------
!               xra = 0.

               xant(i,k) = MAX( xant(i,k) + xra, small_value )
               xnh3(i,k) = MAX( xnh3(i,k) - xra, small_value )
               xhno3(i,k)= MAX( xhno3(i,k)- xra, small_value )
            end if

!-----------------------------------------------------------------
!      ... Washout SO2, SO4 and NH3
!-----------------------------------------------------------------
            xso4(i,k)  = MAX( xso4(i,k), small_value )
            xant(i,k)  = MAX( xant(i,k), small_value )
            xnh3(i,k)  = MAX( xnh3(i,k), small_value )
            xso2(i,k)  = MAX( xso2(i,k), small_value )
         end do
      end do

!==============================================================
!       ... Update the mixing ratios
!==============================================================
      do k = 1,plev
         if( so2_ndx > 0 ) then
            qin(:,k,so2_ndx) =  MAX( xso2(:,k), small_value )
         end if
         if( so4_ndx > 0 ) then
            qin(:,k,so4_ndx) =  MAX( xso4(:,k), small_value )
         end if
         if( h2o2_ndx > 0 ) then
            qin(:,k,h2o2_ndx) =  MAX( xh2o2(:,k), small_value ) 
         end if
         if( nh3_ndx > 0 ) then
            qin(:,k,nh3_ndx) =  MAX( xnh3(:,k), small_value )
         end if
         if( nh4no3_ndx > 0 ) then
            qin(:,k,nh4no3_ndx) =  MAX( xant(:,k), small_value )
         end if
         if( hno3_ndx > 0 ) then
            qin(:,k,hno3_ndx) =  MAX( xhno3(:,k), small_value )
         end if
      end do 
 
      end subroutine SETSOX

      end module MO_SETSOX_MOD
