MODULE polysvp_mod

use fms_mod,             only : FATAL, error_mesg, write_version_number
use sat_vapor_pres_mod,  only : compute_qs, sat_vapor_pres_init
use constants_mod,       only : rdgas, rvgas, hlv, hls, tfreeze, cp_air
use mg_const_mod,        only : mg_const_init, mg_pr
use strat_cloud_utilities_mod,      &
                         only : strat_cloud_utilities_init, &
                                atmos_state_type, cloud_state_type,&
                                strat_nml_type

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

interface polysvp_l
   module procedure polysvp_l1d, polysvp_l3d
end interface polysvp_l

interface polysvp_i
   module procedure polysvp_i1d, polysvp_i3d
end interface polysvp_i

PUBLIC  polysvp_l, polysvp_i, polysvp_init, polysvp_end, &
        compute_qs_a, compute_qs_x1
private compute_qs_x2

!-------------------------------------------------------------------------
!----version number-------------------------------------------------------
Character(len=128) :: Version = '$Id: polysvp.F90,v 20.0 2013/12/13 23:22:03 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'


!---------------------------------------------------------------------
real, parameter :: d622 = rdgas / rvgas
real, parameter :: d378 = 1. - d622

logical         :: module_is_initialized = .false.



CONTAINS



!#########################################################################

SUBROUTINE polysvp_init

      IF (module_is_initialized) return

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number (version, tagname)

!------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!------------------------------------------------------------------------
      CALL sat_vapor_pres_init
      call strat_cloud_utilities_init
      call mg_const_init

!------------------------------------------------------------------------
!    mark this module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .TRUE.


END SUBROUTINE polysvp_init



!#########################################################################

SUBROUTINE compute_qs_a (idim, jdim, kdim, Nml, Atmos_state, Cloud_state)

!-------------------------------------------------------------------------
INTEGER,                INTENT(IN )   :: idim, jdim, kdim
type(strat_nml_type),   intent(in)    :: Nml
type(atmos_state_type), intent(inout) :: Atmos_state
type(cloud_state_type), intent(inout) :: Cloud_state

!-----------------------------------------------------------------------
!----local variables

      INTEGER :: i, j, k

!-------------------------------------------------------------------------
!    compute qs and associated parameters, using different algorithms as
!    appropriate.    
!-------------------------------------------------------------------------
      if (Nml%do_pdf_clouds .and. Nml%pdf_org) then
        if (Nml%super_ice_opt .LT. 1) then

!-------------------------------------------------------------------------
!    if using pdf clouds and not allowing supersaturation:
!-------------------------------------------------------------------------
          call compute_qs (Atmos_state%T_in - ((hlv*Cloud_state%ql_in +  &
                                         hls*Cloud_state%qi_in)/cp_air),  &
                           Atmos_state%pfull, Atmos_state%qs,   &
                           dqsdT=Atmos_state%dqsdT, esat=Atmos_state%esat0)
          Atmos_state%gamma = Atmos_state%dqsdT*(min(1.,  &
            max(0., 0.05*(Atmos_state%T_in - ((hlv*Cloud_state%ql_in +  &
               hls*Cloud_state%qi_in)/cp_air) - tfreeze + 20.)))*hlv +     &
               min(1., max(0., 0.05*(tfreeze - Atmos_state%T_in +   &
                          ((hlv*Cloud_state%ql_in +   &
                              hls*Cloud_state%qi_in)/cp_air))))*hls)/cp_air
        else

!-------------------------------------------------------------------------
!    if using pdf clouds and allowing supersaturation:
!-------------------------------------------------------------------------
          call error_mesg ( 'strat_cloud_mod', &
             'super_ice_opt > 1 and pdf_opt currently not supported', FATAL)
!    this code is pattern for what might be used when this option supported:
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
               !call compute_qs_x1 (T(i,j,k), pfull(i,j,k), qs(i,j,k), &
               !                    qsl(i,j,k), qsi(i,j,k), dqsdT(i,j,k), &
               !                    gamma(i,j,k))
              end do
            end do
          end do 
          Atmos_state%gamma = Atmos_state%dqsdT*(min(1.,  &
               max(0., 0.05*(Atmos_state%T_in - tfreeze + 20.)))*hlv +     &
                       min(1., max(0., 0.05*(tfreeze -   &
                                           Atmos_state%T_in)))*hls)/cp_air
        end if
      else

!-------------------------------------------------------------------------
!    if not using pdf clouds and not allowing supersaturation:
!-------------------------------------------------------------------------
        if (Nml%super_ice_opt .LT. 1) then

!------------------------------------------------------------------------
!    Calculate saturation specific humidity and its temperature 
!    derivative, and thermal conductivity plus vapor diffusivity factor.
!
!    FOR ORIGINAL SCHEME These are calculated according to the formulas:
!
!    (1)  qs   = d622*esat/ [pfull  -  (1.-d622)*esat]
!
!    (2) dqsdT = d622*pfull*(desat/dT)/[pfull-(1.-d622)*esat]**2.
!
!    (3) gamma = (L/cp) * dqsdT
!       
!       where d622 = rdgas/rvgas; esat = saturation vapor pressure;
!       and desat/dT is the temperature derivative of esat.
!       Note that in the calculation of gamma, 
!
!            {             hlv          for T > tfreeze             }
!       L =  { 0.05*(T-tfreeze+20.)*hlv + 0.05*(tfreeze-T)*hls      }
!            {                          for tfreeze-20.< T < tfreeze}
!            {             hls          for T < tfreeze-20.         }
!
!       This linear form is chosen because at tfreeze-20. es = esi, and
!       at tfreeze, es = esl, with linear interpolation in between.
!
!
!    (5) chi = 2.21 E-05 (m*m)/s  * (1.E+05)/pfull
!
!        where p is the pressure in Pascals.
!    
!
!       Note that qs, dqsdT, and gamma do not have their proper values
!       until all of the following code has been executed.  That
!       is qs and dqsdT are used to store intermediary results
!       in forming the full solution.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    calculate water saturated vapor pressure from table and store 
!    temporarily in the variable esat0.
!------------------------------------------------------------------------
          call compute_qs   &
                    (Atmos_state%T_in, Atmos_state%pfull, Atmos_state%qs,&
                           dqsdT=Atmos_state%dqsdT, esat=Atmos_state%esat0)
          Atmos_state%gamma = Atmos_state%dqsdT*(min(1.,   &
            max(0., 0.05*(Atmos_state%T_in - tfreeze + 20.)))*hlv +     &
             min(1., max(0., 0.05*(tfreeze - Atmos_state%T_in)))*hls)/cp_air
        else

!-------------------------------------------------------------------------
!    if not using pdf clouds and allowing supersaturation:
!-------------------------------------------------------------------------
          call compute_qs_x1     &
               (idim, jdim, kdim, Atmos_state%T_in, Atmos_state%pfull, &
                     Atmos_state%qs, Atmos_state%qsl, Atmos_state%qsi,   &
                                      Atmos_state%dqsdT, Atmos_state%gamma)
        end if
      end if

!--------------------------------------------------------------------------
!    Calculate relative humidity outside of convective cloud portion of 
!    grid box (U_ca) under the assumption that the temp is uniform across 
!    the gridbox. Define the maximum portion of the remaining grid-box
!    area which can be cloudy while maintaining the grid-box-mean relative
!    humidity (U01). Here qrat is the ratio of the gridbox RH to that in the
!    environment of the convective cloud. 
!--------------------------------------------------------------------------
      IF (Nml%super_ice_opt .LT. 1) THEN

!--------------------------------------------------------------------------
!    when supersaturation not permitted, upper limit of 100% placed on 
!    external relative humidity.
!--------------------------------------------------------------------------
        where (Atmos_state%qrat .gt. 0.)
          Atmos_state%U_ca = min(max(0.,   &
                 (Atmos_state%qv_in/(Atmos_state%qrat*Atmos_state%qs))), 1.)
        elsewhere
          Atmos_state%U_ca = 0.
        end where
        Atmos_state%U01 = Atmos_state%U_ca
      ELSE

!------------------------------------------------------------------------
!    when supersaturation is permitted, no upper limit placed on external
!    relative humidity. Max cloud fraction  depends on whether clouds are
!    either liquid or ice.
!--------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (Atmos_state%qrat(i,j,k) .gt. 0. ) then
                Atmos_state%U_ca(i,j,k) = max(0.,  &
                  (Atmos_state%qv_in(i,j,k)/(Atmos_state%qrat(i,j,k)*  &
                                                 Atmos_state%qs(i,j,k))))
                if (Atmos_State%T_in(i,j,k) .LT. tfreeze) then
                  Atmos_state%U01(i,j,k) = max(0.,  &
                     (Atmos_state%qv_in(i,j,k)/(Atmos_state%qrat(i,j,k)* &
                                                 Atmos_state%qsi(i,j,k))))
                else
                  Atmos_state%U01(i,j,k) = max(0.,  &
                     (Atmos_state%qv_in(i,j,k)/(Atmos_state%qrat(i,j,k)* &
                                                   Atmos_state%qsl(i,j,k))))
                endif
              else
                Atmos_state%U_ca(i,j,k) = 0.
              endif     
            end do
          end do
        end do
      END IF

!-----------------------------------------------------------------------

END SUBROUTINE compute_qs_a



!#########################################################################

SUBROUTINE  polysvp_end

      module_is_initialized = .FALSE.

END SUBROUTINE  polysvp_end


!#########################################################################

function polysvp_l1d (T)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to liquid  by using 
!    function from Goff and Gatch (1946). polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr) :: T,polysvp_l1d
REAL(kind=mg_pr) :: dum

!-------------------------------------------------------------------------
!    Goff Gatch equation, uncertain below -70 C
!-------------------------------------------------------------------------
      polysvp_l1d = 10._mg_pr**(-7.90298_mg_pr*(373.16_mg_pr/t-1._mg_pr) + &
                    5.02808_mg_pr*log10(373.16_mg_pr/t) - &
                    1.3816e-7_mg_pr*(10._mg_pr**(11.344_mg_pr*(1._mg_pr - &
                                          t/373.16_mg_pr)) - 1._mg_pr) + &
                    8.1328e-3_mg_pr*(10._mg_pr**(-3.49149_mg_pr*   &
                            (373.16_mg_pr/t - 1._mg_pr)) - 1._mg_pr) + &
                                          log10(1013.246_mg_pr))*100._mg_pr 

!-------------------------------------------------------------------------


end function polysvp_l1d



!##########################################################################

function polysvp_i1d (T)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to ice by using 
!    function from Goff and Gatch (1946). Polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr) ::  T,polysvp_i1d
REAL(kind=mg_pr) :: dum

!-----------------------------------------------------------------------
! Goff Gatch equation  for ice (good down to -100 C)
!-----------------------------------------------------------------------
      polysvp_i1d = 10._mg_pr**(-9.09718_mg_pr*(273.16_mg_pr/t -   &
                                            1._mg_pr) - 3.56654_mg_pr* &
                     log10(273.16_mg_pr/t) + 0.876793_mg_pr*  &
                                         (1._mg_pr - t/273.16_mg_pr) + &
                                            log10(6.1071_mg_pr))*100._mg_pr

!-------------------------------------------------------------------------

end function polysvp_i1d


!########################################################################

function polysvp_l3d (T, idim, jdim, kdim)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to liquid  by using 
!    function from Goff and Gatch (1946). polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr), dimension(idim,jdim,kdim) :: T,polysvp_l3d
REAL(kind=mg_pr)                            :: dum
integer                                     :: idim,jdim,kdim

!-------------------------------------------------------------------------
!    Goff Gatch equation, uncertain below -70 C
!-------------------------------------------------------------------------
      polysvp_l3d = 10._mg_pr**(-7.90298_mg_pr*(373.16_mg_pr/t-1._mg_pr) + &
                    5.02808_mg_pr*log10(373.16_mg_pr/t) - &
                    1.3816e-7_mg_pr*(10._mg_pr**(11.344_mg_pr*(1._mg_pr - &
                                          t/373.16_mg_pr)) - 1._mg_pr) + &
                    8.1328e-3_mg_pr*(10._mg_pr**(-3.49149_mg_pr*   &
                            (373.16_mg_pr/t - 1._mg_pr)) - 1._mg_pr) + &
                                          log10(1013.246_mg_pr))*100._mg_pr 

!-------------------------------------------------------------------------


end function polysvp_l3d


!#########################################################################

function polysvp_i3d (T, idim, jdim, kdim)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to ice by using 
!    function from Goff and Gatch (1946). Polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr), dimension(idim,jdim,kdim) :: T, polysvp_i3d
REAL(kind=mg_pr)                            :: dum
integer                                     :: idim, jdim, kdim

!-----------------------------------------------------------------------
! Goff Gatch equation  for ice (good down to -100 C)
!-----------------------------------------------------------------------
      polysvp_i3d = 10._mg_pr**(-9.09718_mg_pr*(273.16_mg_pr/t -   &
                                            1._mg_pr) - 3.56654_mg_pr* &
                     log10(273.16_mg_pr/t) + 0.876793_mg_pr*  &
                                         (1._mg_pr - t/273.16_mg_pr) + &
                                            log10(6.1071_mg_pr))*100._mg_pr

!-------------------------------------------------------------------------
   
end function polysvp_i3d


!#########################################################################

SUBROUTINE compute_qs_x1 (idim, jdim, kdim, ttmp, pfull, qs, qs_l, qs_i,&
                          dqsdT, gamma )

!--------------------------------------------------------------------------
integer,                         intent(in)     :: idim, jdim, kdim
REAL, dimension(idim,jdim,kdim), INTENT(IN )    :: ttmp, pfull
REAL, dimension(idim,jdim,kdim), INTENT(INOUT ) :: qs, dqsdT, gamma,  &
                                                   qs_l, qs_i

!--------------------------------------------------------------------------
!---local variables--------------------------------------------------------
      REAL, dimension (idim,jdim,kdim) :: eslt, esit
      integer                          :: i,j,k

!-------------------------------------------------------------------------
!    compute es with respect to liquid and ice.
!-------------------------------------------------------------------------
      eslt = polysvp_l3d (ttmp, idim, jdim, kdim)
      esit = polysvp_i3d (ttmp, idim, jdim, kdim)

!-------------------------------------------------------------------------
!    make sure ice saturation doesn't exceed water saturation at temps 
!    near freezing.
!-------------------------------------------------------------------------
      where (esit .gt. eslt) esit = eslt

!-------------------------------------------------------------------------
!    compute saturation specific humidity over liquid. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!-------------------------------------------------------------------------
      qs_l = pfull - d378*eslt
      qs_l = max(qs_l, eslt) 
      qs_l = d622*eslt/qs_l

!-------------------------------------------------------------------------
!    compute saturation specific humidity over ice. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!    compute saturation specific humidity over liquid.
!-------------------------------------------------------------------------
      qs_i = pfull - d378*esit
      qs_i = max(qs_i, esit) 
      qs_i = d622*esit/qs_i

!--------------------------------------------------------------------------
!    define the appropriate qs, dqsdT and gamma to use, dependent on ambient
!    temperature.
!--------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            IF (ttmp(i,j,k) .GE. tfreeze - 23.16) THEN
              qs(i,j,k) = qs_l(i,j,k)
            ELSE
              qs(i,j,k) = qs_i(i,j,k)
            END IF
!REV#1
!888
!           IF (ttmp(i,j,k) .GE. 233.15) THEN
!           IF (ttmp(i,j,k) .GE. 233.16) THEN
            IF (ttmp(i,j,k) .GE. tfreeze - 40.) THEN
!END REV#1
              dqsdT(i,j,k) = hlv*qs_l(i,j,k)/(rvgas*ttmp(i,j,k)**2)
              gamma(i,j,k) = dqsdT(i,j,k)*hlv/cp_air
            ELSE
              dqsdT(i,j,k) = hls*qs_i(i,j,k)/(rvgas*ttmp(i,j,k)**2)
              gamma(i,j,k) = dqsdT(i,j,k)*hls/cp_air
            END IF
          end do
        end do
     end do

!------------------------------------------------------------------------

END SUBROUTINE compute_qs_x1


!########################################################################

SUBROUTINE compute_qs_x2 (ttmp, pfull, qs, qs_l, qs_i, dqsdT, gamma,  &
                          ql, qi, qmin )

REAL,    INTENT(IN )    :: ttmp, pfull, qmin, ql, qi
REAL,    INTENT(INOUT ) :: qs,  dqsdT, gamma,  qs_l, qs_i

      REAL :: eslt, esit, ifrac

!-------------------------------------------------------------------------
!    determine ice fraction of total condensate.
!-------------------------------------------------------------------------
      if (qi .gt. qmin) then
        ifrac = qi/(ql + qi)
      else
        ifrac = 0.
      end if

!-------------------------------------------------------------------------
!    compute es with respect to liquid and ice.
!-------------------------------------------------------------------------
      eslt = polysvp_l1d(ttmp)
      esit = polysvp_i1d(ttmp)

!-------------------------------------------------------------------------
!    make sure ice saturation doesn't exceed water saturation at temps 
!    near freezing.
!-------------------------------------------------------------------------
      IF (esit .GT. eslt) esit = eslt 

!-------------------------------------------------------------------------
!    compute saturation specific humidity over liquid. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!-------------------------------------------------------------------------
      qs_l = pfull - d378*eslt
      qs_l = max(qs_l, eslt) 
      qs_l = d622*eslt/qs_l

!-------------------------------------------------------------------------
!    compute saturation specific humidity over ice. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!    compute saturation specific humidity over liquid.
!-------------------------------------------------------------------------
      qs_i = pfull - d378*esit
      qs_i = max(qs_i, esit) 
      qs_i = d622*esit/qs_i

!--------------------------------------------------------------------------
!    define the appropriate qs, dqsdT and gamma to use, dependent on ambient
!    temperature and ice fraction.
!--------------------------------------------------------------------------
      IF (ttmp .GE. tfreeze - 23.16 .and.  ifrac .LT. 0.9 ) THEN
        qs =  qs_l
      ELSE
        qs =  qs_i
      END IF
!REV#2
!888
!     IF (ttmp .GE. 233.15  .and.  ifrac .LT. 0.9 ) THEN
!     IF (ttmp .GE. 233.16  .and.  ifrac .LT. 0.9 ) THEN
      IF (ttmp .GE. tfreeze - 40.  .and.  ifrac .LT. 0.9 ) THEN
!END REV#2
        dqsdT =  hlv*qs_l/(rvgas*ttmp**2)
        gamma = dqsdT * hlv /cp_air
      ELSE
        dqsdT =  hls*qs_i/(rvgas*ttmp**2)
        gamma = dqsdT * hls/ cp_air
      END IF

!-------------------------------------------------------------------------

END SUBROUTINE compute_qs_x2


!#########################################################################

END MODULE polysvp_mod
