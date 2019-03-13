! ----------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of MOM.                                                                 
!                                                                      
! MOM is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
!
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Surface fCO2 calculation
!</OVERVIEW>
!
!<DESCRIPTION>
! Calculate the fugacity of CO2 at the surface in thermodynamic
! equilibrium with the current alkalinity (Alk) and total dissolved
! inorganic carbon (DIC) at a particular temperature and salinity
! using an initial guess for the total hydrogen
! ion concentration (htotal)
!</DESCRIPTION>
!

module ocmip2_co2calc_mod

implicit none

private

public  :: ocmip2_co2calc
public  :: ocmip2_co2_alpha

character(len=128) :: version = '$Id: ocmip2_co2calc.F90,v 20.0 2013/12/14 00:09:44 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

contains


!#######################################################################
! <SUBROUTINE NAME="ocmip2_co2_alpha">
!
! <DESCRIPTION>
!       Calculate CO2 solubility, alpha, from
! temperature (t) and salinity (s).
!
! INPUT
!
!       isd        = first i-limit of the arrays with halo
!       ied        = last i-limit of the arrays with halo
!       jsd        = first j-limit of the arrays with halo
!       jed        = last j-limit of the arrays with halo
!       isc        = first i-limit of the arrays for computation
!       iec        = last i-limit of the arrays for computation
!       jsc        = first j-limit of the arrays for computation
!       jec        = last j-limit of the arrays for computation
!
!       t          = temperature (degrees C)
!
!       s          = salinity (PSU)
!
!       mask       = land mask array (0.0 = land)
!
! OUTPUT
!       alpha      = Solubility of CO2 for air? (mol/kg/atm unless scaled)
!
! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!
!     - Models may carry tracers in mol/m^3 (on a per volume basis)
!
!     - Conversely, this routine, which was written by observationalists
!       (C. Sabine and R. Key), passes input arguments in umol/kg  
!       (i.e., on a per mass basis)
!
!     - Thus, if the input or output units need to be changed from/to mol/m^3
!       then set scale to an appropriate value. For example, if the model
!       uses mol/m^3, then scale should be set to something like 1.0/1024.5
!       to convert from mol/m^3 to mol/kg.
!
! </DESCRIPTION>
subroutine ocmip2_co2_alpha(isd, jsd, isc, iec, jsc, jec, t, s, mask, alpha, scale)

integer, intent(in)                     :: isd
integer, intent(in)                     :: jsd
integer, intent(in)                     :: isc
integer, intent(in)                     :: iec
integer, intent(in)                     :: jsc
integer, intent(in)                     :: jec
real, dimension(isd:,jsd:), intent(in)  :: t
real, dimension(isd:,jsd:), intent(in)  :: s
real, dimension(isd:,jsd:), intent(in)  :: mask
real, dimension(isc:,jsc:), intent(out) :: alpha
real, intent(in), optional              :: scale

integer :: i, j
real    :: log100
real    :: tk
real    :: tk100
real    :: tk1002
real    :: logtk
real    :: ff
real    :: scale_factor

! set the scale factor for unit conversion
if (present(scale)) then
  scale_factor = scale
else
  scale_factor = 1.0
endif

! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Once calculated, the constants are stored and passed in the common
! block "const". The original version of this code was based on
! the code by Dickson in Version 2 of "Handbook of Methods for the
! Analysis of the Various Parameters of the Carbon Dioxide System
! in Seawater", DOE, 1994 (SOP No. 3, p25-26).
!
! Derive simple terms used more than once
log100 = log(100.0)
do j = jsc, jec
  do i = isc, iec

    if (mask(i,j) .ne. 0.0) then

      tk     = 273.15 + t(i,j)
      tk100  = tk / 100.0
      tk1002 = tk100 * tk100
      logtk  = log(tk)

      ! f = k0(1-pH2O)*correction term for non-ideality
      
      ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6
      ! values)
      ff = exp(-162.8301 + 218.2968 / tk100  +                                  &
               90.9241 * (logtk - log100) - 1.47696 * tk1002 +                  &
               s(i,j) * (0.025695 - 0.025225 * tk100 + 0.0049867 * tk1002))

      !      Should we be using K0 or ff for the solubility here?
      !         convert to output units
      alpha(i,j) = ff / scale_factor
    else
      alpha(i,j) = 0.0
    endif
  enddo
enddo

return

end subroutine  ocmip2_co2_alpha
! </SUBROUTINE> NAME="ocmip2_co2_alpha"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_co2calc">
!
! <DESCRIPTION>
!       Calculate co2* from total alkalinity and total CO2 at
! temperature (t) and salinity (s).
! It is assumed that init_ocmip2_co2calc has already been called with
! the T and S to calculate the various coefficients.
!
! INPUT
!
!       isd        = first i-limit of the arrays with halo
!       ied        = last i-limit of the arrays with halo
!       jsd        = first j-limit of the arrays with halo
!       jed        = last j-limit of the arrays with halo
!       isc        = first i-limit of the arrays for computation
!       iec        = last i-limit of the arrays for computation
!       jsc        = first j-limit of the arrays for computation
!       jec        = last j-limit of the arrays for computation
!
!       mask       = land mask array (0.0 = land)
!
!       t          = temperature (degrees C)
!
!       s          = salinity (PSU)
!
!       dic_in     = total inorganic carbon (mol/kg unless scaled) 
!                    where 1 T = 1 metric ton = 1000 kg
!
!       ta_in      = total alkalinity (eq/kg unless scaled) 
!
!       pt_in      = inorganic phosphate (mol/kg unless scaled) 
!
!       sit_in     = inorganic silicate (mol/kg unless scaled) 
!
!       htotallo   = factor to set lower limit of htotal range
!
!       htotalhi   = factor to set upper limit of htotal range
!
!       htotal     = H+ concentraion
!
! OUTPUT
!       co2star    = CO2*water (kg/kg unless scaled)
!       alpha      = Solubility of CO2 for air? (kg/kg/atm unless scaled)
!       pco2surf   = oceanic pCO2 (ppmv)
!
!       k1         = activity factors for carbonate species
!
!       k2              (see below)
!
!       invtk      = 1/(t+273.15)
!
! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!
!     - Models may carry tracers in mol/m^3 (on a per volume basis)
!
!     - Conversely, this routine, which was written by observationalists
!       (C. Sabine and R. Key), passes input arguments in umol/kg  
!       (i.e., on a per mass basis)
!
!     - Thus, if the input or output units need to be changed from/to mol/m^3
!       then set scale to an appropriate value. For example, if the model
!       uses mol/m^3, then scale should be set to something like 1.0/1024.5
!       to convert from mol/m^3 to mol/kg.
!
! </DESCRIPTION>
subroutine ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec, mask, t, s,   &
     dic_in, ta_in, pt_in, sit_in, htotallo, htotalhi, htotal,                  &
     co2star, co3_ion, alpha, pCO2surf, k1_out, k2_out, invtk_out, scale)

real, parameter :: permeg = 1.e-6
real, parameter :: xacc = 1.0e-10

integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
real, dimension(isd:,jsd:), intent(in)                  :: mask
real, dimension(isd:,jsd:), intent(in)                  :: t
real, dimension(isd:,jsd:), intent(in)                  :: s
real, dimension(isd:,jsd:), intent(in)                  :: dic_in
real, dimension(isd:,jsd:), intent(in)                  :: ta_in
real, dimension(isd:,jsd:), intent(in)                  :: pt_in
real, dimension(isd:,jsd:), intent(in)                  :: sit_in
real, dimension(isc:,jsc:), intent(in)                  :: htotallo
real, dimension(isc:,jsc:), intent(in)                  :: htotalhi
real, dimension(isc:,jsc:), intent(inout)               :: htotal
real, dimension(isc:,jsc:), intent(out), optional       :: co2star
real, dimension(isc:,jsc:), intent(out), optional       :: co3_ion
real, dimension(isc:,jsc:), intent(out), optional       :: alpha
real, dimension(isc:,jsc:), intent(out), optional       :: pCO2surf
real, dimension(isc:,jsc:), intent(out), optional       :: invtk_out
real, dimension(isc:,jsc:), intent(out), optional       :: k1_out
real, dimension(isc:,jsc:), intent(out), optional       :: k2_out
real, intent(in), optional                              :: scale

integer :: i
integer :: j
real    :: log100
real    :: pt
real    :: sit
real    :: ta
real    :: dic
real    :: tk
real    :: invtk
real    :: tk100
real    :: tk1002
real    :: logtk
real    :: is
real    :: is2
real    :: sqrtis
real    :: s2
real    :: sqrts
real    :: s15
real    :: scl
real    :: logf_of_s
real    :: ff
!real    :: k0
real    :: k1
real    :: k2
real    :: kb
real    :: k1p
real    :: k2p
real    :: k3p
real    :: ksi
real    :: kw
real    :: ks
real    :: kf
real    :: bt
real    :: st
real    :: ft
real    :: htotal2
real    :: co2star_internal
real    :: scale_factor

!       set the scale factor for unit conversion
if (present(scale)) then
  scale_factor = scale
else
  scale_factor = 1.0
endif

! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Once calculated, the constants are stored and passed in the common
! block "const". The original version of this code was based on
! the code by Dickson in Version 2 of "Handbook of Methods for the
! Analysis of the Various Parameters of the Carbon Dioxide System
! in Seawater", DOE, 1994 (SOP No. 3, p25-26).
!
! Derive simple terms used more than once
log100 = log(100.0)
do j = jsc, jec
  do i = isc, iec

    if (mask(i,j) .ne. 0.0) then

      tk        = 273.15 + t(i,j)
      tk100     = tk / 100.0
      tk1002    = tk100 * tk100
      invtk     = 1.0 / tk
      logtk    = log(tk)
      is        = 19.924 * s(i,j) / (1000.0 - 1.005 * s(i,j))
      is2       = is * is
      sqrtis    = sqrt(is)
      s2        = s(i,j) * s(i,j)
      sqrts     = sqrt(s(i,j))
      s15       = sqrts ** 3
      scl       = s(i,j) / 1.80655
      logf_of_s = log(1.0 - 0.001005 * s(i,j))

      ! k0 from Weiss 1974

      !k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +                 &
                !s(i,j) * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))


      ! k1 = [H][HCO3]/[H2CO3]
      ! k2 = [H][CO3]/[HCO3]

      ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
      k1 = 10.0 ** (-(3670.7 * invtk - 62.008 + 9.7944 * logtk -                &
                    0.0118 * s(i,j) + 0.000116 * s2))

      k2 = 10.0 ** (-(1394.7 * invtk + 4.777 - 0.0184 * s(i,j) + 0.000118 * s2))

      ! kb = [H][BO2]/[HBO2]

      ! Millero p.669 (1995) using data from Dickson (1990)
      kb = exp((-8966.90 - 2890.53 * sqrts - 77.942 * s(i,j) +                  &
                1.728 * s15 - 0.0996 * s2) * invtk +                            &
               (148.0248 + 137.1942 * sqrts + 1.62142 * s(i,j)) +               &
               (-24.4344 - 25.085 * sqrts - 0.2474 * s(i,j)) * logtk +          &
               0.053105 * sqrts * tk)

      ! k1p = [H][H2PO4]/[H3PO4]

      ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
      k1p = exp(-4576.752 * invtk + 115.525 - 18.453 * logtk +                  &
                (-106.736 * invtk + 0.69171) * sqrts +                          &
                (-0.65643 * invtk - 0.01844) * s(i,j))

      ! k2p = [H][HPO4]/[H2PO4]

      ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
      k2p = exp(-8814.715 * invtk + 172.0883 - 27.927 * logtk +                 &
                (-160.340 * invtk + 1.3566) * sqrts +                           &
                (0.37335 * invtk - 0.05778) * s(i,j))

      ! k3p = [H][PO4]/[HPO4]

      ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
      k3p = exp(-3070.75 * invtk - 18.141 +                                     &
                (17.27039 * invtk + 2.81197) * sqrts +                          &
                (-44.99486 * invtk - 0.09984) * s(i,j))

      ! ksi = [H][SiO(OH)3]/[Si(OH)4]

      ! Millero p.671 (1995) using data from Yao and Millero (1995)
      ksi = exp(-8904.2 * invtk + 117.385 - 19.334 * logtk +                    &
                (-458.79 * invtk + 3.5913) * sqrtis +                           &
                (188.74 * invtk - 1.5998) * is +                                &
                (-12.1652 * invtk + 0.07871) * is2 + logf_of_s)

      ! kw = [H][OH]
      
      ! Millero p.670 (1995) using composite data
      kw = exp(-13847.26 * invtk + 148.9652 - 23.6521 * logtk +                 &
               (118.67 * invtk - 5.977 + 1.0495 * logtk) * sqrts - 0.01615 * s(i,j))

      ! ks = [H][SO4]/[HSO4]

      ! Dickson (1990, J. chem. Thermodynamics 22, 113)
      ks = exp(-4276.1 * invtk + 141.328 - 23.093 * logtk +                     &
               (-13856.0 * invtk + 324.57 - 47.986 * logtk) * sqrtis +          &
               (35474.0 * invtk - 771.54 + 114.723 * logtk) * is -              &
               2698.0 * invtk * sqrtis ** 3 + 1776.0 * invtk * is2 + logf_of_s)

      ! kf = [H][F]/[HF]
      
      ! Dickson and Riley (1979) -- change pH scale to total
      kf = exp(1590.2 * invtk - 12.641 + 1.525 * sqrtis + logf_of_s +           &
               log(1.0 + (0.1400 / 96.062) * scl / ks))

      ! Calculate concentrations for borate, sulfate, and fluoride
      
      ! Uppstrom (1974)
      bt = 0.000232 / 10.811 * scl

      ! Morris & Riley (1966)
      st = 0.14 / 96.062 * scl

      ! Riley (1965)
      ft = 0.000067 / 18.9984 * scl

      ! set some stuff to pass back, if requested
      if (present(k1_out)) then
        k1_out(i,j) = k1
      endif

      if (present(k2_out)) then
        k2_out(i,j) = k2
      endif

      if (present(invtk_out)) then
        invtk_out(i,j)  = invtk
      endif

      ! Possibly convert input in mol/m^3 -> mol/kg 
      sit = sit_in(i,j) * scale_factor
      ta  = ta_in(i,j)  * scale_factor
      dic = dic_in(i,j) * scale_factor
      pt  = pt_in(i,j)  * scale_factor

!
!***********************************************************************
!
! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.
!
! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate htotal. In this case we use the
! Newton-Raphson "safe" method taken from "Numerical Recipes"
! (function "rtsafe.f" with error trapping removed).
!
! As currently set, this procedure iterates about 12 times. The x1
! and x2 values set below will accomodate ANY oceanographic values.
! If an initial guess of the pH is known, then the number of
! iterations can be reduced to about 5 by narrowing the gap between
! x1 and x2. It is recommended that the first few time steps be run
! with x1 and x2 set as below. After that, set x1 and x2 to the
! previous value of the pH +/- ~0.5. The current setting of xacc will
! result in co2star accurate to 3 significant figures (xx.y). Making
! xacc bigger will result in faster convergence also, but this is not
! recommended (xacc of 10**-9 drops precision to 2 significant
! figures).
!

      htotal(i,j) = drtsafe(htotalhi(i,j)*htotal(i,j),  &
           htotallo(i,j)*htotal(i,j),                   &
           dic, ta, pt, sit, k1, k2, k1p, k2p, k3p,     &
           bt, ft, st, kb, kw, kf, ks, ksi, xacc)

      ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
      ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
      ! Convert units of output arguments
      !      Note: co2star is calculated in
      !            mol/kg within this routine 
      !
      htotal2 = htotal(i,j) * htotal(i,j)
      co2star_internal = dic * htotal2 / (htotal2 + k1 * (htotal(i,j) + k2)) / scale_factor
      if (present(co2star)) then
        co2star(i,j) = co2star_internal
      endif
      if (present(co3_ion)) then
        co3_ion(i,j) = co2star_internal * k1 * k2 / htotal2
      endif
      !ph = -log10(htotal(i,j))

      ! f = k0(1-pH2O)*correction term for non-ideality

      ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6
      !                values)
      if (present(alpha) .or. present(pco2surf)) then

        ff = exp(-162.8301 + 218.2968 / tk100  +                                        &
                 90.9241 * (logtk - log100) - 1.47696 * tk1002 +                        &
                 s(i,j) * (0.025695 - 0.025225 * tk100 + 0.0049867 * tk1002))

        ! Add two output arguments for storing pCO2surf
        ! Should we be using K0 or ff for the solubility here?
        ! Convert pCO2surf from atm to uatm
        ! possibly convert output units
        if (present(alpha)) then
          alpha(i,j) = ff / scale_factor
        endif
        if (present(pco2surf)) then
          pCO2surf(i,j) = co2star_internal / (ff / scale_factor) / permeg
        endif

      endif

    else

      if (present(co2star)) then
        co2star(i,j) = 0.0
      endif
      if (present(k1_out)) then
        k1_out(i,j) = 0.0
      endif
      if (present(k2_out)) then
        k2_out(i,j) =0.0
      endif
      if (present(invtk_out)) then
        invtk_out(i,j)  = 0.0
      endif
      if (present(alpha)) then
        alpha(i,j) = 0.0
      endif
      if (present(pco2surf)) then
        pCO2surf(i,j) = 0.0
      endif

    endif

  enddo
enddo

return

end subroutine  ocmip2_co2calc
! </SUBROUTINE> NAME="ocmip2_co2calc"


!#######################################################################
! <FUNCTION NAME="drtsafe">
!
! <DESCRIPTION>
!       File taken from Numerical Recipes. Modified  R. M. Key 4/94
! </DESCRIPTION>

function drtsafe(x1, x2, dic, ta, pt, sit, k1, k2, k1p, k2p, k3p,       &
                 bt, ft, st, kb, kw, kf, ks, ksi, xacc)

real                    :: drtsafe
real, intent(in)        :: x1
real, intent(in)        :: x2
real, intent(in)        :: dic
real, intent(in)        :: ta
real, intent(in)        :: pt
real, intent(in)        :: sit
real, intent(in)        :: k1
real, intent(in)        :: k2
real, intent(in)        :: k1p
real, intent(in)        :: k2p
real, intent(in)        :: k3p
real, intent(in)        :: bt
real, intent(in)        :: ft
real, intent(in)        :: st
real, intent(in)        :: kb
real, intent(in)        :: kw
real, intent(in)        :: kf
real, intent(in)        :: ks
real, intent(in)        :: ksi
real, intent(in)        :: xacc

integer, parameter      :: maxit = 100

integer :: j
real    :: fl
real    :: df
real    :: fh
real    :: swap
real    :: xl
real    :: xh
real    :: dxold
real    :: dx
real    :: f
real    :: temp

call ta_iter_1(x1, fl, df, dic, ta, pt, sit, k1, k2,    &
     k1p, k2p, k3p, bt, ft, st, kb, kw, kf, ks, ksi)
call ta_iter_1(x2, fh, df, dic, ta, pt, sit, k1, k2,    &
     k1p, k2p, k3p, bt, ft, st, kb, kw, kf, ks, ksi)
if(fl .lt. 0.0) then
  xl=x1
  xh=x2
else
  xh=x1
  xl=x2
  swap=fl
  fl=fh
  fh=swap
end if
drtsafe=0.5*(x1+x2)
dxold=abs(x2-x1)
dx=dxold
call ta_iter_1(drtsafe, f, df, dic, ta, pt, sit, k1, k2,        &
     k1p, k2p, k3p, bt, ft, st, kb, kw, kf, ks, ksi)
do j=1,maxit
  if (((drtsafe-xh)*df-f)*((drtsafe-xl)*df-f) .ge. 0.0 .or.     &
      abs(2.0*f) .gt. abs(dxold*df)) then
    dxold=dx
    dx=0.5*(xh-xl)
    drtsafe=xl+dx
    if (xl .eq. drtsafe) then
      return
    endif
  else
    dxold=dx
    dx=f/df
    temp=drtsafe
    drtsafe=drtsafe-dx
    if (temp .eq. drtsafe) then
      return
    endif
  end if
  if (abs(dx) .lt. xacc) then
    return
  endif
  call ta_iter_1(drtsafe, f, df, dic, ta, pt, sit, k1, k2,      &
       k1p, k2p, k3p, bt, ft, st, kb, kw, kf, ks, ksi)
  if(f .lt. 0.0) then
    xl=drtsafe
    fl=f
  else
    xh=drtsafe
    fh=f
  end if
enddo

return

end function drtsafe
! </FUNCTION> NAME="drtsafe"


!#######################################################################
! <SUBROUTINE NAME="ta_iter_1">
!
! <DESCRIPTION>
! This routine expresses TA as a function of DIC, htotal and constants.
! It also calculates the derivative of this function with respect to 
! htotal. It is used in the iterative solution for htotal. In the call
! "x" is the input value for htotal, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhtotal
! </DESCRIPTION>

subroutine ta_iter_1(x, fn, df, dic, ta, pt, sit, k1, k2,       &
     k1p, k2p, k3p, bt, ft, st, kb, kw, kf, ks, ksi)

real, intent(in)        :: x
real, intent(out)       :: fn
real, intent(out)       :: df
real, intent(in)        :: dic
real, intent(in)        :: ta
real, intent(in)        :: pt
real, intent(in)        :: sit
real, intent(in)        :: k1
real, intent(in)        :: k2
real, intent(in)        :: k1p
real, intent(in)        :: k2p
real, intent(in)        :: k3p
real, intent(in)        :: bt
real, intent(in)        :: ft
real, intent(in)        :: st
real, intent(in)        :: kb
real, intent(in)        :: kw
real, intent(in)        :: kf
real, intent(in)        :: ks
real, intent(in)        :: ksi

real    :: x2
real    :: x3
real    :: k12
real    :: k12p
real    :: k123p
real    :: c
real    :: a
real    :: a2
real    :: da
real    :: b
real    :: b2
real    :: db

x2 = x*x
x3 = x2*x
k12 = k1*k2
k12p = k1p*k2p
k123p = k12p*k3p
c = 1.0 + st/ks
a = x3 + k1p*x2 + k12p*x + k123p
a2 = a*a
da = 3.0*x2 + 2.0*k1p*x + k12p
b = x2 + k1*x + k12
b2 = b*b
db = 2.0*x + k1

!     fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
fn = k1*x*dic/b + 2.0*dic*k12/b + bt/(1.0 + x/kb) + kw/x +              &
     pt*k12p*x/a + 2.0*pt*k123p/a + sit/(1.0 + x/ksi) -                 &
     x/c - st/ (1.0 + ks/x/c) - ft/(1.0 + kf/x) -                       &
     pt*x3/a - ta

!     df = dfn/dx
df = ((k1*dic*b) - k1*x*dic*db)/b2 - 2.0*dic*k12*db/b2 -                &
     bt/kb/(1.0+x/kb)**2 - kw/x2 + (pt*k12p*(a - x*da))/a2 -            &
     2.0*pt*k123p*da/a2 - sit/ksi/(1.0+x/ksi)**2 - 1.0/c +              &
     st*(1.0 + ks/x/c)**(-2)*(ks/c/x2) + ft*(1.0 + kf/x)**(-2)*kf/x2 -  &
     pt*x2*(3.0*a-x*da)/a2

return

end subroutine  ta_iter_1
! </SUBROUTINE> NAME="ta_iter_1"

end module  ocmip2_co2calc_mod
