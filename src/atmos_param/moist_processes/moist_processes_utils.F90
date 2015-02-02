
                    module moist_proc_utils_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             relaxed arakawa-schubert
!             donner deep convection
!             large-scale condensation
!             stratiform prognostic cloud scheme 
!             rel humidity cloud scheme 
!             diagnostic cloud scheme 
!             lin cloud microphysics
!             betts-miller convective adjustment
!
!-----------------------------------------------------------------------

use  sat_vapor_pres_mod, only: compute_qs, lookup_es
use    time_manager_mod, only: time_type
use    diag_manager_mod, only: send_data
use       constants_mod, only: RDGAS, RVGAS

implicit none
private

!------------------ private and public data/interfaces -----------------

public   capecalcnew, tempavg, column_diag, rh_calc

private  column_diag_1, column_diag_2, column_diag_3

interface column_diag
  module procedure column_diag_1, column_diag_2, column_diag_3
end interface column_diag


!--------------------- version number ----------------------------------
character(len=128) :: &
version = '$Id: moist_processes_utils.F90,v 19.0 2012/01/06 20:10:44 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------

real, public, allocatable, dimension(:,:,:) :: pmass

                             contains

!#######################################################################

      subroutine tempavg (pdepth,phalf,temp,tsnow,mask)

!-----------------------------------------------------------------------
!
!    computes a mean atmospheric temperature for the bottom
!    "pdepth" pascals of the atmosphere.
!
!   input:  pdepth     atmospheric layer in pa.
!           phalf      pressure at model layer interfaces
!           temp       temperature at model layers
!           mask       data mask at model layers (0.0 or 1.0)
!
!   output:  tsnow     mean model temperature in the lowest
!                      "pdepth" pascals of the atmosphere
!
!-----------------------------------------------------------------------
      real, intent(in)  :: pdepth
      real, intent(in) , dimension(:,:,:) :: phalf,temp
      real, intent(out), dimension(:,:)   :: tsnow
      real, intent(in) , dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
 real, dimension(size(temp,1),size(temp,2)) :: prsum, done, pdel, pdep
 real  sumdone
 integer  k
!-----------------------------------------------------------------------

      tsnow=0.0; prsum=0.0; done=1.0; pdep=pdepth

      do k=size(temp,3),1,-1

         if (present(mask)) then
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*mask(:,:,k)*done(:,:)
         else
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*done(:,:)
         endif

         where ((prsum(:,:)+pdel(:,:))  >  pdep(:,:))
            pdel(:,:)=pdepth-prsum(:,:)
            done(:,:)=0.0
            pdep(:,:)=0.0
         endwhere

         tsnow(:,:)=tsnow(:,:)+pdel(:,:)*temp(:,:,k)
         prsum(:,:)=prsum(:,:)+pdel(:,:)

         sumdone=sum(done(:,:))
         if (sumdone < 1.e-4) exit

      enddo

         tsnow(:,:)=tsnow(:,:)/prsum(:,:)

!-----------------------------------------------------------------------

      end subroutine tempavg

!#######################################################################
!cape calculation.
                                                                                
subroutine capecalcnew(kx,p,phalf,cp,rdgas,rvgas,hlv,kappa,tin,rin,&
                                avgbl,cape,cin)
                                                                                
!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp          specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron,
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!    avgbl       if true, the parcel is averaged in theta and r up to its LCL
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set
!                to zero)
!
!    Algorithm:
!    Start with surface parcel.
!    Calculate the lifting condensation level (uses an analytic formula and a
!       lookup table).
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.
      implicit none
      integer, intent(in)                    :: kx
      logical, intent(in)                    :: avgbl
      real, intent(in), dimension(:)         :: p, phalf, tin, rin
      real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp
      real, intent(out)                      :: cape, cin
                                                                                
      integer            :: k, klcl, klfc, klzb
      logical            :: nocape
      real, dimension(kx)   :: tp, rp
      real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
                               a, b, dtdlnp, &
                               plcl, plzb

      pstar = 1.e5
                                                                                
      nocape = .true.
      cape = 0.
      cin = 0.
      plcl = 0.
      plzb = 0.
      klfc = 0
      klcl = 0
      klzb = 0
      tp(1:kx) = tin(1:kx)
      rp(1:kx) = rin(1:kx)
                                                                                
! start with surface parcel
      t0 = tin(kx)
      r0 = rin(kx)
! calculate the lifting condensation level by the following:
! are you saturated to begin with?
      call lookup_es(t0,es)
      rs = rdgas/rvgas*es/p(kx)
      if (r0.ge.rs) then
! if you're already saturated, set lcl to be the surface value.
         plcl = p(kx)
! the first level where you're completely saturated.
         klcl = kx
! saturate out to get the parcel temp and humidity at this level
! first order (in delta T) accurate expression for change in temp
         tp(kx) = t0 + (r0 - rs)/(cp/hlv + hlv*rs/rvgas/t0**2.)
         call lookup_es(tp(kx),es)
         rp(kx) = rdgas/rvgas*es/p(kx)
      else
! if not saturated to begin with, use the analytic expression to calculate the
! exact pressure and temperature where you?re saturated.
         theta0 = tin(kx)*(pstar/p(kx))**kappa
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! The right hand side of this is only a function of temperature, therefore
! this is put into a lookup table to solve for temperature.
         if (r0.gt.0.) then
            value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas) 
            call lcltabl(value,tlcl)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
! just in case plcl is very high up
            if (plcl.lt.p(1)) then
               plcl = p(1)
               tlcl = theta0*(plcl/pstar)**kappa
               write (*,*) 'hi lcl'
            end if
            k = kx
         else
! if the parcel sp hum is zero or negative, set lcl to 2nd to top level
            plcl = p(2)
            tlcl = theta0*(plcl/pstar)**kappa
!            write (*,*) 'zero r0', r0
            do k=2,kx
               tp(k) = theta0*(p(k)/pstar)**kappa
               rp(k) = 0.
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            end do
            go to 11
         end if
! calculate the parcel temperature (adiabatic ascent) below the LCL.
! the mixing ratio stays the same
         do while (p(k).gt.plcl)
            tp(k) = theta0*(p(k)/pstar)**kappa
            call lookup_es(tp(k),es)
            rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            k = k-1
         end do
! first level where you're saturated at the level
         klcl = k
         if (klcl.eq.1) klcl = 2
! do a saturated ascent to get the parcel temp at the LCL.
! use your 2nd order equation up to the pressure above.
! moist adaibat derivatives: (use the lcl values for temp, humid, and
! pressure)
         a = kappa*tlcl + hlv/cp*r0
         b = hlv**2.*r0/cp/rvgas/tlcl**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
! second order in p (RK2)
! first get temp halfway up
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
         a = kappa*tp(klcl) + hlv/cp*rp(klcl)
         b = hlv**2./cp/rvgas*rp(klcl)/tp(klcl)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
!         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
!                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp/ &
!                   (1. + b)
! second order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
!             p(klcl)/plcl))**2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/p(klcl)
!         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
! CAPE/CIN stuff
         if ((tp(klcl).lt.tin(klcl)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(klcl) - &
                 tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(klcl) - &
                  tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = klcl
            endif
         end if
      end if
! then, start at the LCL, and do moist adiabatic ascent by the first order
! scheme -- 2nd order as well
      do k=klcl-1,1,-1
         a = kappa*tp(k+1) + hlv/cp*rp(k+1)
         b = hlv**2./cp/rvgas*rp(k+1)/tp(k+1)**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
! second order in p (RK2)
! first get temp halfway up
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
         a = kappa*tp(k) + hlv/cp*rp(k)
         b = hlv**2./cp/rvgas*rp(k)/tp(k)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
!         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - &
!               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp*rp(k+1)/(1. + b)
! second order in p

!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
!             p(k)/p(k+1)))**2.
! if you're below the lookup table value, just presume that there's no way
! you could have cape and call it quits
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/p(k)
         if ((tp(k).lt.tin(k)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
         elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
! if you have CAPE, and it's your first time being negatively buoyant,
! then set the level of zero buoyancy to k+1, and stop the moist ascent
            klzb = k+1
            go to 11
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = k
            endif
         end if
      end do
 11   if(nocape) then
! this is if you made it through without having a LZB
! set LZB to be the top level.
         plzb = p(1)
         klzb = 0
         klfc = 0
         cin = 0.
         tp(1:kx) = tin(1:kx)
         rp(1:kx) = rin(1:kx)
      end if
!      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
!      write (*,*) 'tp, rp new', tp, rp
!       write (*,*) 'tp, new', tp
!       write (*,*) 'tin new', tin
!       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
      end subroutine capecalcnew


!#######################################################################
! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.
!
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! the RHS is tabulated for the control amount of moisture, hence the 
! division by es00 on the LHS

! Gives the values of the temperature for the following range:
!   starts with -23, is uniformly distributed up to -10.4.  There are a
! total of 127 values, and the increment is .1.
!
      implicit none
      real, intent(in)     :: value
      real, intent(out)    :: tlcl
      integer              :: ival
      real, dimension(127) :: lcltable
      real                 :: v1, v2
                                                                                
      data lcltable/   1.7364512e+02,   1.7427449e+02,   1.7490874e+02, &
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02, &
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02, &
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02, &
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02, &
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02, &
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02, &
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02, &
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02, &
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02, &
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02, &
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02, &
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02, &
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02, &
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02, &
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02, &
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02, &
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02, &
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02, &
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02, &
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02, &
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02, &
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02, &
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02, &
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02, &
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02, &
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02, &
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02, &
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02, &
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02, &
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02, &
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02 /
                                                                                
      v1 = value
      if (value.lt.-23.0) v1 = -23.0
      if (value.gt.-10.4) v1 = -10.4
      ival = floor(10.*(v1 + 23.0))
      v2 = -230. + ival
      v1 = 10.*v1
      tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(ival+2)
                                                                                
      end subroutine lcltabl


!#######################################################################
subroutine column_diag_1 (id_diag, is, js, Time, val1, c_val1, temp_in)
  integer, intent(in)                       :: id_diag, is, js
  real, intent(in)                          :: c_val1
  type(time_type), intent(in)               :: Time
  real, dimension(:,:,:), intent(in)        :: val1
  real, dimension(:,:), optional, intent(inout) :: temp_in
!local
  real, dimension(size(val1,1), size(val1,2)) :: temp
  integer :: k
  logical :: used
  integer :: ie, je

  if (present(temp_in)) then
    temp = temp_in
  else
    temp = 0.
  endif

  ie = is + size(val1,1) -1
  je = js + size(val1,2) -1

  do k = 1,size(val1,3)
    temp(:,:) = temp(:,:) + c_val1 * val1(:,:,k) * pmass(is:ie,js:je,k)
  enddo
  used = send_data (id_diag, temp, Time, is, js)

  if (present(temp_in)) temp_in=temp

end subroutine column_diag_1


!#######################################################################
subroutine column_diag_2 (id_diag, is, js, Time, val1, c_val1, val2, c_val2, temp_in)
  integer, intent(in)                       :: id_diag, is, js
  real, intent(in)                          :: c_val1, c_val2
  type(time_type), intent(in)               :: Time
  real, dimension(:,:,:), intent(in)        :: val1, val2
  real, dimension(:,:), optional, intent(inout) :: temp_in
!local
  real, dimension(size(val1,1), size(val1,2)) :: temp
  integer :: k
  logical :: used
  integer :: ie, je

  if (present(temp_in)) then
    temp = temp_in
  else
    temp = 0.
  endif

  ie = is + size(val1,1) -1
  je = js + size(val1,2) -1
  do k = 1,size(val1,3)
    temp(:,:) = temp(:,:) + (c_val1 * val1(:,:,k) + &
                             c_val2 * val2(:,:,k) ) * pmass(is:ie,js:je,k)
  enddo
  used = send_data (id_diag, temp, Time, is, js)

  if (present(temp_in)) temp_in=temp

end subroutine column_diag_2


!#######################################################################
subroutine column_diag_3 (id_diag, is, js, Time, val1, c_val1, val2, c_val2, val3, c_val3, temp_in)
  integer, intent(in)                       :: id_diag, is, js
  real, intent(in)                          :: c_val1, c_val2, c_val3
  type(time_type), intent(in)               :: Time
  real, dimension(:,:,:), intent(in)        :: val1, val2, val3
  real, dimension(:,:), optional, intent(inout) :: temp_in
!local
  real, dimension(size(val1,1), size(val1,2)) :: temp
  integer :: k
  logical :: used
  integer :: ie, je

  if (present(temp_in)) then
    temp = temp_in
  else
    temp = 0.
  endif

  ie = is + size(val1,1) -1
  je = js + size(val1,2) -1
  do k = 1,size(val1,3)
    temp(:,:) = temp(:,:) + (c_val1 * val1(:,:,k) + &
                             c_val2 * val2(:,:,k) + &
                             c_val3 * val3(:,:,k) ) * pmass(is:ie,js:je,k)
  enddo
  used = send_data (id_diag, temp, Time, is, js)

  if (present(temp_in)) temp_in=temp

end subroutine column_diag_3

  
!#######################################################################
subroutine rh_calc(pfull,T,qv,RH,do_simple,MASK, do_cmip)
  
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity. 
!       This is calculated according to the formula:
! 
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!     
!       Where epsilon = RDGAS/RVGAS = d622
!     
!       and where 1- epsilon = d378
!     
!       Note that RH does not have its proper value 
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        IMPLICIT NONE
        
        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK
        logical, intent(in), optional :: do_cmip
        logical, intent(in) :: do_simple
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat
      
        real, parameter :: d622 = RDGAS/RVGAS
        real, parameter :: d378 = 1.-d622

! because Betts-Miller uses a simplified scheme for calculating the relative humidity
        if (do_simple) then
          call lookup_es(T, esat)
          RH(:,:,:) = pfull(:,:,:)
          RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:))  !limit where pfull ~ esat
          RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))
        else
          if (present(do_cmip)) then
            call compute_qs (T, pfull, rh, q=qv,  &
                                          es_over_liq_and_ice = .true.)
             RH(:,:,:)=qv(:,:,:)/RH(:,:,:)
          else
            call compute_qs (T, pfull, rh, q=qv)
            RH(:,:,:)=qv(:,:,:)/RH(:,:,:)
          endif
        endif

        !IF MASK is present set RH to zero
        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)

END SUBROUTINE rh_calc


!#######################################################################

                 end module moist_proc_utils_mod

