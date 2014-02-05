 
!VERSION NUMBER:
!  $Name: tikal $
!  $Id: wet_deposition_0D.F90,v 19.0 2012/01/06 20:09:12 fms Exp $

!<SUBROUTINE NAME = "wet_deposition_0D">
!<TEMPLATE>
!CALL wet_deposition_0D( Henry_constant, Henry_variable, &
!                        frac_in_cloud, alpha_r, alpha_s, &
!                        T, p0, p1, rho_air, &
!                        cloud, precip, &
!                        tracer, Lgas, Laerosol, Lice, &
!                        delta_tracer )
!</TEMPLATE>
subroutine wet_deposition_0D( Henry_constant, Henry_variable, &
                              frac_in_cloud, alpha_r, alpha_s, &
                              T, p0, p1, rho_air, &
                              cloud, rain, snow, &
                              tracer, Lgas, Laerosol, Lice, &
                              delta_tracer )
implicit none
!      
!<OVERVIEW>
! Routine to calculate the fraction of tracer removed by wet deposition
!</OVERVIEW>
!
!<IN NAME="T" TYPE="real">
!   Temperature (K)
!</IN>
!<IN NAME="p0" TYPE="real">
!   Pressure (Pa) at layer closer to surface
!</IN>
!<IN NAME="p1" TYPE="real">
!   Pressure (Pa) at layer farther from surface
!</IN>
!<IN NAME="rho_air" TYPE="real">
!   Air density (kg/m3)
!</IN>
!<IN NAME="cloud" TYPE="real">
!   Cloud amount (liquid+ice) (kg/kg)
!</IN>
!<IN NAME="rain" TYPE="real">
!   Precipitation increment (rain) (kg/m3)
!</IN>
!<IN NAME="snow" TYPE="real">
!   Precipitation increment (snow) (kg/m3)
!</IN>
!<IN NAME="tracer" TYPE="real">
!   The tracer field (tracer units)
!</IN>
!<IN NAME="Lgas" TYPE="logical">
!   Is tracer a gas?
!</IN>
!<IN NAME="Laerosol" TYPE="logical">
!   Is tracer an aerosol?
!</IN>
!<IN NAME="Lice" TYPE="logical">
!   Is tracer removed by snow (or only by rain)?
!</IN>
!<OUT NAME="delta_tracer" TYPE="real">
!   The change (increment) of the tracer field due to wet deposition (tracer units)
!/OUT>
!<DESCRIPTION>
! Schemes allowed here are:
!
! 1) Removal according to Henry's Law. This law states that the ratio of the concentation in 
!    cloud water and the partial pressure in the interstitial air is a constant. In this 
!    instance, the units for Henry's constant are kg/L/Pa (normally it is M/L/Pa)
!    Parameters for a large number of species can be found at
!    http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!
! 2) Aerosol removal, using specified in-cloud tracer fraction

! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "wet_deposition","henry","henry=XXX, dependence=YYY"
! "wet_deposition","henry_below","henry=XXX, dependence=YYY"
!     where XXX is the Henry's constant for the tracer in question
!       and YYY is the temperature dependence of the Henry's Law constant.
!
! "wet_deposition","aerosol","frac_incloud=XXX"
! "wet_deposition","aerosol_below","frac_incloud=XXX"
!     where XXX is the in-cloud fraction of the aerosol tracer
!</PRE>

!</DESCRIPTION>

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
real,             intent(in)                     :: Henry_constant, Henry_variable, &
                                                    frac_in_cloud, alpha_r, alpha_s
real,             intent(in)                     :: T, p0, p1, rho_air
real,             intent(in)                     :: cloud, rain, snow
real,             intent(in)                     :: tracer
logical,          intent(in)                     :: Lgas, Laerosol, Lice
real,             intent(out)                    :: delta_tracer

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
real :: &
      Htemp, xliq, n_air, pwt, pmid, precip
real :: &
      temp_factor, scav_factor, &
      w_h2o, beta, f_a, in_temp
real, parameter :: &
      GRAV = 9.80,              &  ! acceleration due to gravity [m/s2]
      RDGAS = 287.04,           &  ! gas constant for dry air [J/kg/deg]
      AVOGNO = 6.023000E+23,    &  ! Avogadro's number
      inv298p15 = 1./298.15,    &  ! 1/K
      cm3_2_m3 = 1.e-6             ! m3/cm3
real, parameter :: mw_air = 28.96440E-03 ! molar mass of air (kg/mole)
real, parameter :: mw_h2o = 18.0E-03  ! molar mass of H2O (kg/mole)

!-----------------------------------------------------------------------

delta_tracer = 0.

pmid = 0.5 * (p0+p1)         ! Pa
pwt     = ( p0 - p1 )/GRAV   ! kg/m2

if( Lgas .or. Laerosol ) then
!++lwh
! units = VMR
!
! Henry_constant (mole/L/Pa) = [X](aq) / Px(g) 
! where [X](aq) is the concentration of tracer X in precipitation (mole/L)
!       Px(g) is the partial pressure of the tracer in the air (Pa)
!
! VMR (total) = VMR (gas) + VMR (aq)
!             = VMR (gas) + [X] * L
!
! where L = cloud liquid amount (kg H2O/mole air)
!
! Using Henry's Law, [X] = H * Px = H * VMR(gas) * Pfull
!
! So, VMR (total) =  VMR(gas) * [ 1 + H * Pfull * L ]
! 
! VMR(gas) = VMR(total) / [1 + H * Pfull * L]
!
! [X] = H * Pfull * VMR(total) / [ 1 + H * Pfull * L]
!
! Following Giorgi and Chameides, JGR, 90(D5), 1985, the first-order loss
! rate constant (s^-1) of X due to wet deposition equals:
!
! k = W_X / n_X
!
! where W_x = the loss rate (molec/cm3/s), and n_X = the number density (molec/cm3)
! 
! W_X = [X] * W_H2O / (55 mole/L)
! n_x = VMR(total) * n_air (molec/cm3) = VMR(total) * P/(kT) * 1E-6 m3/cm3
! 
! where P = atmospheric pressure (Pa)
!       k = Boltzmann's constant = 1.38E-23 J/K
!       T = temperature (K)
!       W_H2O = removal increment of water (molec/cm3)
! 
!             [X] * W_H2O / 55         
! So, k = ------------------------------
!         VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O    H * VMR(total) * P / [ 1 + H * P *L ]
!       = ----- * ---------------------------------------
!          55          VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O     H * kT * 1E6
!       = ----- *  -------------    
!          55      1 + H * P * L 
!
!         W_H2O     1     1     H * P * L
!       = ----- * ----- * - * -------------
!          55     n_air   L   1 + H * P * L
!
! where W_H2O = precip (kg/m3) * (AVOGNO/mw_h2o) (molec/kg) * 1E-6 m3/cm3
!
   if( (Lgas .and. Henry_constant > 0.) .or. Laerosol ) then
      if (Lice) then
         precip = rain+snow
      else
         precip = rain
      end if
      in_temp = 0.

      scav_factor = 0.0
      xliq = MAX( cloud * mw_air, 0. ) ! (kg H2O)/(mole air)
      n_air = rho_air * (AVOGNO/mw_air) * cm3_2_m3 ! molec/cm3
      if (Lgas) then
! Calculate the temperature dependent Henry's Law constant
         temp_factor = 1/T-inv298p15
         Htemp = Henry_constant * exp( Henry_variable*temp_factor )
         f_a = Htemp * pmid * xliq
         scav_factor = f_a / ( 1.+f_a )
      else if (Laerosol) then
         scav_factor = frac_in_cloud
      end if
      if (precip > 0. .and. xliq > 0.) then
         w_h2o = precip * (AVOGNO/mw_h2o) * cm3_2_m3 ! molec/cm3
         beta = w_h2o * mw_h2o  / (n_air * xliq)   ! fraction of condensed water removed
         beta = MAX(MIN(beta,1.),0.)
         in_temp = beta * scav_factor              ! fraction of tracer removed
      end if

!     wdep_in = - in_temp*tracer*pwt
!     dt_temp = 1. - exp( -in_temp*dt ) ! fractional loss/timestep
!     tracer_dt = dt_temp / dt !+ve loss frequency (1/sec)
      delta_tracer = in_temp ! fraction of tracer removed
!--lwh
   endif 

end if

! Now multiply by the tracer mixing ratio to get the actual tendency.
! tracer_dt = MIN( MAX(tracer_dt, 0.0E+00), 0.5/dt)
if (tracer > 0.) then
   delta_tracer = delta_tracer*tracer
else
   delta_tracer = 0.
end if

! Output diagnostics in kg/m2/s (if MMR) or mole/m2/s (if VMR)
! if(trim(units) .eq. 'mmr') then
!    diag_scale = 1.
! else if(trim(units) .eq. 'vmr') then
!    diag_scale = mw_air ! kg/mole
! else
!    write(*,*) ' Tracer number =',n,' tracer_name=',tracer_name
!    write(*,*) ' scheme=',text_in_scheme
!    write(*,*) ' control=',control
!    write(*,*) ' scheme=',scheme
!    write(*,*) 'Please check field table'
!    write(*,*) 'tracers units =',trim(units),'it should be either  mmr or vmr!'
!  <ERROR MSG="Unsupported tracer units" STATUS="FATAL">
!     Tracer units must be either VMR or MMR
!  </ERROR>
!    call error_mesg('wet_deposition', 'Unsupported tracer units.', FATAL )
! end if

! if(trim(cloud_param) == 'donner') then
!    if (id_tracer_wdep_donin(n) > 0 ) then
!        used = send_data ( id_tracer_wdep_donin(n), wdep_in/diag_scale, Time, is_in=is, js_in=js)
!    endif
!    if(id_tracer_wdep_donin_dt(n) > 0) then
!       used = send_data ( id_tracer_wdep_donin_dt(n), in_temp, Time, is_in=is, js_in=js, ks_in=1)
!    endif
!    if(id_tracer_wdep_don_dt(n) > 0) then
!       used = send_data ( id_tracer_wdep_don_dt(n), dt_temp/dt, Time, is_in=is, js_in=js, ks_in=1)
!    endif
! endif

end subroutine wet_deposition_0D
!</SUBROUTINE>
