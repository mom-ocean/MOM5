MODULE ice_nucl_mod

!-------------------------------------------------------------------------
! based on Liu et al., J Clim, 2007, p4526 ff, Liu and Penner, 2005
! see also Salzmann et al., 2010, to be submitted to ACPD
! and cjg's aer_ccn_nucl_wpdf 
! compute IN activation by integrating over assumed subgrid-scale PDF of w
!-------------------------------------------------------------------------

use mpp_mod,           only : input_nml_file
use fms_mod,           only : error_mesg, FATAL, mpp_pe, mpp_root_pe, &
                              open_namelist_file, check_nml_error, &
                              close_file, write_version_number, &
                              file_exist, stdlog
use aer_ccn_act_k_mod, only : ghquad, dlocate
use aerosol_params_mod,only : aerosol_params_init, rho_sulf, sigma_sulf,  &
                              rho_bc, sigma_bc, Nfact_du1, Nfact_du2,  &
                              Nfact_du3, Nfact_du4, Nfact_du5, Nfact_ss1, &
                              Nfact_ss2, Nfact_ss3, Nfact_ss4, Nfact_ss5

implicit none
private

!-------------------------------------------------------------------------
!--interfaces-------------------------------------------------------------

public ice_nucl_wpdf, ice_nucl_wpdf_init, ice_nucl_wpdf_end
private ice_nucl_k, fast, slow, bc_het

!------------------------------------------------------------------------
!----version number------------------------------------------------------
Character(len=128) :: Version = '$Id: ice_nucl.F90,v 20.0 2013/12/13 23:24:20 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!------------------------------------------------------------------------
!--namelist--------------------------------------------------------------

integer     :: dust_opt = 1             ! option for dust particle 
                                        ! activation expression to be used
logical     :: do_het = .true.          ! allow heterogeneous nucleation ?
logical     :: use_dust_instead_of_bc =   &
                                .false. ! use dust instead of black
                                        ! carbon for heterogeneous 
                                        ! nucleation ?
logical     :: limit_immersion_frz =      &
                                .false. ! limit balck carbon heterogeneous
                                        ! nucleation to number of droplets 
                                        ! present  ??
logical     :: limit_rhil = .false.     !
logical     :: do_ice_nucl_ss_wpdf =      &
                                .false. ! use seasalt particles for 
                                        ! homogeneous nucleation ?
! ---> h1g
logical     :: retain_ice_nucl_bug = .true.
! <--- h1g
!------------------------------------------------------------------------
! note that in-situ nucleation at cold temperatures most likely results
! in small sulfate (homogeneous nucleation). Heterogeneous nucleation, 
! on the other hand, might require larger particles, why it probably
! makes sense to assume that soot has a larger mean diameter.  
!------------------------------------------------------------------------
real        :: d_sulf = 0.04e-6         ! mean diameter of sulfate particles
                                        !  for homogeneous nucleation
!alternate values:
!              d_sulf = 0.1e-6          ! Haywood and Ramaswamy, JGR, 1998
!              d_sulf = 0.04e-6         ! Barahona and Nenes
!              d_sulf = 0.02e-6

real        :: d_bc = 0.0236e-6         ! mean diameter of black carbon 
                                        ! particles for heterogeneous 
                                        ! nucleation
!alternate values:
!              d_bc = 0.0236e-6         ! Hess et al., bams 98
!              d_bc = 0.07e-6           !Pueschel et al., GRL,1992 
!              d_bc = 0.04e-6 

real        :: rh_crit_het = 1.2        !
real        :: dust_surf = 0.5          ! 
real        :: dust_frac_min = 0.0      ! 
real        :: dust_frac_max = 1.e15    !
real        :: dust_frac =1.0           ! fraction of total available dust 
                                        ! to use for nucleation
!RSH: Should this be 150. or 1.5 (150 %) ????
real        :: rh_dust_max = 150.       ! maximum ice supersaturation in the
                                        ! presence of dust


namelist / ice_nucl_nml /  dust_opt, do_het, use_dust_instead_of_bc, &
                           limit_immersion_frz, limit_rhil, &
                           do_ice_nucl_ss_wpdf, d_sulf, &
                           d_bc, rh_crit_het, dust_surf, dust_frac_min, &
                           dust_frac_max, dust_frac, rh_dust_max, retain_ice_nucl_bug

integer, parameter :: npoints = 64     ! # for Gauss-Hermite quadrature
real, parameter    :: wp2_eps = 0.0001 ! w variance threshold
real, parameter    :: wmin =  0.0      ! min w for ccn_nucl
real, parameter    :: wmax = 10.0      ! max w for ccn_nucl
REAL, PARAMETER    :: pi = 3.14
REAL, PARAMETER    :: naer_min = 1.0 
REAL, PARAMETER    :: naer_max = 1.e11
real(kind=8),       dimension(npoints) :: x, w
real, parameter    :: w_het_thresh = 0.001 ! min vert vel for which 
                                           ! heterogeneous nucleation occurs
real, parameter    :: w_hom_thresh = 0.001 ! min vert vel for which 
                                           ! homogeneous nucleation occurs


logical            :: module_is_initialized = .false.


!-----------------------------------------------------------------------

contains

!#########################################################################

SUBROUTINE ice_nucl_wpdf_init

      integer :: unit, io, ierr, logunit
     
!------------------------------------------------------------------------
      IF (module_is_initialized) return

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=ice_nucl_nml, iostat=io)
      ierr = check_nml_error(io,'ice_nucl_nml')
#else
      if ( file_exist('input.nml')) then
 
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=ice_nucl_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'ice_nucl_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------- write version and namelist to standard log ------------
 
        call write_version_number(version, tagname)
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
                       write ( logunit, nml=ice_nucl_nml )
 

!------------------------------------------------------------------------
!    be sure other needed modules are initialized.
!------------------------------------------------------------------------
      call aerosol_params_init

!------------------------------------------------------------------------
!    initialize arrays with abscissas and weights for integration.
!------------------------------------------------------------------------
      call ghquad( npoints, x, w )

      module_is_initialized = .TRUE.

!------------------------------------------------------------------------

END SUBROUTINE  ice_nucl_wpdf_init


!########################################################################

SUBROUTINE ice_nucl_wpdf (T, rhi, rhl, wm, wp2, zfull, totalmass, imass, &
                          tym, tyms, crystal, drop, hom, rh_crit_1d,  &
                          rh_crit_min_1d, ni_sulf, ni_dust, ni_bc )

!------------------------------------------------------------------------
!    subroutine ice_nucl_wpdf computes ice nucleation assuming a normal 
!    distribution of w given by its mean (wm) and second moment (wp2)
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
integer, intent(in)    ::  tym, tyms
real,    intent(in)    ::  zfull,  totalmass(tym), imass(tyms), T, wm, wp2,&
                           rhl, rhi, drop
real,    intent(out)   ::  rh_crit_1d, rh_crit_min_1d, ni_sulf, ni_dust,  &
                           ni_bc, crystal, hom

!--------------------------------------------------------------------------
!---local variables

      logical lintegrate
      integer ia, ib
      real wtmp
      real(kind=8) :: tmp, a, b, sum1, sum2, sum1a, sum1x1, sum1x2, sum1x3
      integer iw

!-------------------------------------------------------------------------
!    initialize variable used to catch minimum.
!-------------------------------------------------------------------------
      rh_crit_min_1d = 1.e9

!-------------------------------------------------------------------------
!    determine whether integration is needed to compute number of nucleated
!    crystals. lintegrate = .true. indicates that numerical integration is 
!    to be performed.
!-------------------------------------------------------------------------
      if (wp2 .gt. wp2_eps) then

!-------------------------------------------------------------------------
!    integration bounds: from wmin to wmax (0 to 10 m/s)
!-------------------------------------------------------------------------
        tmp = 1.0d0/sqrt(2.0*wp2) 
        a = (wmin - wm)*tmp
        b = (wmax - wm)*tmp

!-------------------------------------------------------------------------
!    locate indices within integration bounds
!-------------------------------------------------------------------------
        call dlocate( x, npoints, a, ia )
        call dlocate( x, npoints, b, ib )

!-------------------------------------------------------------------------
!    ia (ib) is zero if a (b) is smaller than the lowest abscissa.
!    in that case, start the integration with the first abscissa.
!-------------------------------------------------------------------------
        ia = ia +1 
        ia = min(max(ia,1),size(x))
        ib = min(max(ib,1),size(x))
        if (ib .gt. ia) then
          lintegrate = .true.
        else
          lintegrate = .false.
        endif
      else
        lintegrate = .false.
      endif

!-------------------------------------------------------------------------
!    compute number of nucleated crystals.
!-------------------------------------------------------------------------
      if (lintegrate ) then

!-------------------------------------------------------------------------
!    perform integration
!-------------------------------------------------------------------------
        if (T - 273.15  .le. -35. ) then

!-------------------------------------------------------------------------
!    integration at temps below -35
!-------------------------------------------------------------------------
          sum1 = 0.0d0
          sum1a = 0.0d0
          sum2 = 0.0d0
          sum1x1 = 0.0d0
          sum1x2 = 0.0d0
          sum1x3 = 0.0d0
          tmp = sqrt(2.0*wp2)
          do iw=ia,ib
            wtmp = tmp * x(iw) + wm
            call ice_nucl_k (zfull,  T, rhi, rhl, wtmp, totalmass(1),      &
                             imass, tyms, crystal,  drop, hom, rh_crit_1d, &
                             ni_sulf, ni_dust, ni_bc ) 
            sum1 = sum1 + w(iw)*crystal
            sum1x1 = sum1x1 + w(iw) * ni_sulf
            sum1x2 = sum1x2 + w(iw) * ni_dust
            sum1x3 = sum1x3 + w(iw) * ni_bc
            sum1a = sum1a +  w(iw)*rh_crit_1d
            sum2 = sum2 + w(iw)
            rh_crit_min_1d = MIN(rh_crit_min_1d, rh_crit_1d)
          enddo

!-------------------------------------------------------------------------
!    normalize over the distribution.
!-------------------------------------------------------------------------
          crystal = sum1 / sum2
          rh_crit_1d = sum1a / sum2
          ni_sulf = sum1x1 / sum2
          ni_dust = sum1x2 / sum2
          ni_bc = sum1x3 / sum2

!------------------------------------------------------------------------
!    integration at temps above -35. results here are not a function of
!    the vertical velocity pdf, so only one call need be made to ice_nucl_k.
!-------------------------------------------------------------------------
        else
          call ice_nucl_k (zfull, T, rhi, rhl, wm, totalmass(1),       &
                           imass, tyms, crystal,  drop, hom, rh_crit_1d, &
                           ni_sulf, ni_dust, ni_bc ) 

          rh_crit_min_1d = MIN(rh_crit_min_1d, rh_crit_1d)
        endif

!-------------------------------------------------------------------------
!    no integration, use single point evaluation
!-------------------------------------------------------------------------
      else
        call ice_nucl_k (zfull, T, rhi, rhl, wm, totalmass(1), imass, &
                         tyms, crystal,  drop, hom, rh_crit_1d, ni_sulf,  &
                         ni_dust, ni_bc  )
        rh_crit_min_1d = rh_crit_1d
      endif

!------------------------------------------------------------------------



END SUBROUTINE ice_nucl_wpdf



!#########################################################################

SUBROUTINE  ice_nucl_wpdf_end

      module_is_initialized = .FALSE.

END SUBROUTINE  ice_nucl_wpdf_end




!#########################################################################

SUBROUTINE ice_nucl_k (zfull, T1, rhi_in, rhl_in, W1, TotalMass,      &
                       imass, tyms, Ni, drop, hom, rh_crit_1d,  &
                       ni_sulf, ni_dust, ni_bc) 

!-------------------------------------------------------------------------
integer,              intent(in)    :: tyms          
real,                 intent(in)    :: TotalMass
real, dimension(tyms),intent(in)    :: imass
real,                 intent(in)    :: zfull, T1, rhi_in, W1, rhl_in, drop
real,                 intent(inout) :: hom
real,                 intent(out)   :: Ni, rh_crit_1d, ni_sulf,  &
                                       ni_dust, ni_bc

!--------------------------------------------------------------------------
!--local variables-----

      REAL    :: tc, rhl_thresh, A, B, C, nsulf, nss, naer, nbc, ndu,  &
                 ndu_l, nbccrit,  rhi, rhl, rhid, Sat_max
      LOGICAL :: do_hom

!-----------------------------------------------------------------------
!    place an upper limit on the input relative humidities with respect to
!    liquid and ice.
!-----------------------------------------------------------------------
      rhi = MIN(rhi_in, 2.)
      rhl = MIN(rhl_in, 2.)

!-------------------------------------------------------------------------
!    initialize output fields.
!-------------------------------------------------------------------------
      Ni_bc = 0.
      Ni_sulf = 0.
      Ni_dust = 0.
      rh_crit_1d = rh_crit_het

!-------------------------------------------------------------------------
!    define celsius temperature.
!-------------------------------------------------------------------------
      tc = T1 - 273.15

!------------------------------------------------------------------------
!    calculate ice nucleation at temps below -35C. different procedures are
!    invoked at temps above and below -35C. 
!------------------------------------------------------------------------
      IF ( tc .LE. -35. ) THEN  

!--------------------------------------------------------------------------
!    if heterogeneous nucleation is requested and the vertical velocity
!    exceeds a threshold:
!--------------------------------------------------------------------------
        IF (do_het .and. w1 .GE. w_het_thresh) THEN 

!--------------------------------------------------------------------------
!    use either dust or black carbon as a nucleator, dependent on the nml
!    logical variable use_dust_instead_of_bc.
!--------------------------------------------------------------------------
          IF ( use_dust_instead_of_bc) THEN

!-------------------------------------------------------------------------
!    dust is activated.  units of #/cm^3. Use only dust_frac of the total
!    dust for nucleation.
!-------------------------------------------------------------------------
            ndu =  1.e-6*(Nfact_du1*imass(8) + Nfact_du2*imass(9) + &
                          Nfact_du3*imass(10) + Nfact_du4*imass(11) + &
                          Nfact_du5*imass(12))
            nbc = dust_frac * ndu
          ELSE

!-------------------------------------------------------------------------
!    black carbon is activated. 1.e-6 is conversion from from m^-3 to cm^-3.
!-------------------------------------------------------------------------
            nbc = MIN(imass(6)*1.e-6*6./(rho_bc*pi*d_bc**3)* &
                       exp(-9./2. * (log(sigma_bc))**2), 1.e10)
          ENDIF

!-------------------------------------------------------------------------
!    define the relative humidity threshold for heterogeneous nucleation. 
!    (rhl_thresh). if the rh over liquid exceeds this value, define the 
!    critical number of black carbon particles needed to produce 
!    heterogeneous nucleation (nbccrit). if both thresholds are exceeded,
!    then set a flag indicating homogeneous nucleation is not active, and
!    call bc_het to calculate the number of activated ice nuclei.
!-------------------------------------------------------------------------
          rhl_thresh = 0.0073 * Tc**2 + 1.466 * Tc + 131.74 ! 4.4 L+P05 
          if (1.e2*rhl .GE. rhl_thresh) THEN  ! only heterogeneous nucl.
            nbccrit = exp ( ( 12.884 * log(w1) - 67.69 - Tc ) /  &
                                           (1.4938* log (w1)  + 10.41 ) )
            IF (nbc .GT. nbccrit) THEN ! only heterogeneous nucl.
              do_hom = .FALSE.
              CALL bc_het(Ni_bc, nbc, w1, tc)

!-------------------------------------------------------------------------
!    *rhi/rhl bec. rhi_thresh to be returned to calling program.
!-------------------------------------------------------------------------
              rh_crit_1d = 1.e-2*rhl_thresh*rhi_in/rhl_in

!-------------------------------------------------------------------------
!    limit by  max. # of droplets, if so requested in nml.
!-------------------------------------------------------------------------
              IF (limit_immersion_frz) THEN
                Ni_bc = MIN(Ni_bc, drop)
              END IF

!-------------------------------------------------------------------------
!    convert to  m-3
!-------------------------------------------------------------------------
              Ni_bc = MAX(1.e6*Ni_bc, 0.)

!-------------------------------------------------------------------------
!    if insufficient black carbon particles for heterogenous nucleation, 
!    set flag to allow homogeneous nucleation.
!-------------------------------------------------------------------------
            else
              do_hom = .TRUE.
            endif  

!-------------------------------------------------------------------------
!    if insufficient humidty for heterogenous nucleation, 
!    set flag to allow homogeneous nucleation.
!-------------------------------------------------------------------------
          else
            do_hom = .TRUE.
          endif  

!-------------------------------------------------------------------------
!    if heterogeneous nucleation was not requested, or the vertical velocity
!    is inadequate for heterogenous nucleation, set flag to allow 
!    homogeneous nucleation.
!-------------------------------------------------------------------------
        else
          do_hom = .TRUE.
        endif 

!-------------------------------------------------------------------------
!     calculate homogeneous nucleation.
!-------------------------------------------------------------------------
        IF (do_hom .and. w1 .GE. w_hom_thresh) THEN  

!------------------------------------------------------------------------
!    calculate relative humidity threshold for homogeneous nucleation.
!------------------------------------------------------------------------
          A = 6.e-4 * LOG(W1) + 6.6e-3
! ---> h1g, 2012-06-29, B=6.e-2 * LOG(W1) + 1.052 from 
!   (1) Liu, X., and J. E. Penner, 2005: Ice nucleation parameterization
!                for global models. Meteor. Z., 14, 499-514. (2005)
!   (2) Liu, Xiaohong, Joyce E. Penner, Steven J. Ghan, Minghuai Wang, 2007: 
!                Inclusion of Ice Microphysics in the NCAR Community Atmospheric 
!                Model Version 3 (CAM3). J. Climate, 20, 4526-4547. doi: http://dx.doi.org/10.1175/JCLI4264.1 
!   (3) M. Salzmann1,*, Y. Ming2, J.-C. Golaz2, P. A. Ginoux2, H. Morrison3, A. Gettelman3, M. Kr¨amer4, 
!                and L. J. Donner, Two-moment bulk stratiform cloud microphysics 
!                in the GFDL AM3 GCM: description, evaluation, and sensitivity tests, ACP (2010)
!   change from
!         B = 6.e-3 * LOG(W1) + 1.052
!          to
        if( retain_ice_nucl_bug ) then
          B = 6.e-3 * LOG(W1) + 1.052 
        else
          B = 6.e-2 * LOG(W1) + 1.052
        endif 

! <--- h1g, 2012-06-29

          C = 1.68 * LOG(W1) + 129.35
          rhl_thresh = A * Tc**2 + B * Tc + C 

!------------------------------------------------------------------------
!    *rhi/rhl bec. rhi_thresh to be returned to calling program.
!-----------------------------------------------------------------------
          rh_crit_1d = 1.e-2*rhl_thresh*rhi_in/rhl_in

!-------------------------------------------------------------------------
!    if relative humidity exceeds threshold, calculate the available
!    number of sulfate particles for nucleation.
!    log-normal distr, V=1/6 Pi Dbar^3 N exp(9/2 ln^2 sigma)
!    nsulf in #/cm^3
!-------------------------------------------------------------------------
          IF (1.e2*rhl .GE. rhl_thresh) THEN  
            nsulf = MIN(TotalMass*1.e-6*1.0e-9*1.0e12*6./   &
                             (rho_sulf*pi*d_sulf**3 ) *  &
                                 exp(-9./2. * (log(sigma_sulf))**2), 1.e20)
            naer = nsulf

!-------------------------------------------------------------------------
!    if seasalt particles area desired as ice nuclei, calculate their
!    availability. place upper and lower limits on the total number of 
!    particles available for homogeneous nucleation. 
!-------------------------------------------------------------------------
            IF ( do_ice_nucl_ss_wpdf ) THEN
              nss = 1.e-6*(Nfact_ss1*imass(1) + Nfact_ss2*imass(2) +  &
                           Nfact_ss3*imass(3) + Nfact_ss4*imass(4) +  &
                           Nfact_ss5*imass(5))
              naer = naer + nss
            ENDIF
            naer = MIN(MAX(naer, naer_min), naer_max)

!------------------------------------------------------------------------
!    calculate the number of activated sulfate / seasalt particles. the
!    formula used to determine the number activated is dependent on 
!    temperature and vertical velocity.
!------------------------------------------------------------------------
            IF (tc .GE. 6.07*LOG(W1) - 55.) THEN
              CALL fast (W1, naer, Ni_sulf, tc)
            ELSE
              CALL slow (W1, naer, Ni_sulf, tc)
            END IF 

!------------------------------------------------------------------------
!    convert the activated number to  m-3.
!------------------------------------------------------------------------
            Ni_sulf = MAX(1.e6 * Ni_sulf, 0.)
          END IF 
        END IF 

!-------------------------------------------------------------------------
!    calculate activated ice nuclei when temperature is between -5C and
!    -35C. set the homogeneoyus nucleation flag to .false.
!-------------------------------------------------------------------------
      ELSE IF (tc .GT.  -35. .AND. tc .LT. -5.) THEN
        do_hom = .false.

!-------------------------------------------------------------------------
!    in the presence of dust, ice supersat greater than a certain threshold
!    (rh_dust_max) will not occur.
!-------------------------------------------------------------------------
        rhid = MIN(rhi, rh_dust_max) 

!-------------------------------------------------------------------------
!    optional dust particle activation parameterizations:
!-------------------------------------------------------------------------
        IF (dust_opt .EQ. 1) THEN

!-------------------------------------------------------------------------
!    Yi's factor
!-------------------------------------------------------------------------
          Ni_dust = MIN(MAX(imass(7)/dust_surf, dust_frac_min),   &
                     dust_frac_max)*EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))

!-------------------------------------------------------------------------
!    Meyer's curve.
!-------------------------------------------------------------------------
        ELSE IF (dust_opt .EQ. 2) THEN 
          Ni_dust  =  EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))

!-------------------------------------------------------------------------
!    profile similar to Liu et al. based on !NH! of Minikin et al.
!-------------------------------------------------------------------------
        ELSE IF (dust_opt .EQ. 3) THEN 
          Ni_dust = EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))
          IF (zfull .GT. 1000. .AND. zfull .LE. 7000.) THEN
            Ni_dust  = Ni_dust - 0.9 * Ni_dust/6000.*( zfull - 1000.)
          ELSE IF (zfull .GT. 7000.) THEN
            Ni_dust  = 0.1*Ni_dust
          END IF

!-------------------------------------------------------------------------
!    with scaling factor from Phillips et al.,2008
!-------------------------------------------------------------------------
        ELSE IF (dust_opt .EQ. 4) THEN 
          Ni_dust  = 0.154* MIN(MAX(imass(7)/dust_surf, dust_frac_min),  &
                     dust_frac_max)*EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))
    
!-------------------------------------------------------------------------
!    limit dust nuclei by total amount present. ndu_l = 1.e-3*ndust.
!-------------------------------------------------------------------------
        ELSE IF   ( dust_opt .EQ. 5 ) THEN 
          ndu_l = 1.e-3*(Nfact_du1*imass(8) + Nfact_du2*imass(9)  +  &
                         Nfact_du3*imass(10) + Nfact_du4*imass(11) +  &
                         Nfact_du5*imass(12))
          Ni_dust = MIN(MAX(imass(7)/dust_surf, dust_frac_min),   &
                     dust_frac_max)*EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))
          Ni_dust = MIN (Ni_dust,  ndu_l)

        ELSE IF   ( dust_opt .EQ. 6 ) THEN
          Ni_dust  = 1.5e-10* EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))

! --->h1g, 2012-06-30
! calculate maximum super-saturation Sat_max following 
!    (1) Liu, X., and J. E. Penner, 2005: Ice nucleation parameterization
!                for global models. Meteor. Z., 14, 499-514. (2005)
        ELSE IF   ( dust_opt .EQ. 7 ) THEN
          IF ( use_dust_instead_of_bc) THEN

!-------------------------------------------------------------------------
!    dust is activated.  units of #/cm^3. Use only dust_frac of the total
!    dust for nucleation.
!-------------------------------------------------------------------------
            ndu =  1.e-6*(Nfact_du1*imass(8) + Nfact_du2*imass(9) + &
                          Nfact_du3*imass(10) + Nfact_du4*imass(11) + &
                          Nfact_du5*imass(12))
            nbc = dust_frac * ndu
          ELSE
!-------------------------------------------------------------------------
!    black carbon is activated. 1.e-6 is conversion from from m^-3 to cm^-3.
!-------------------------------------------------------------------------
            nbc =  MIN(imass(6)*1.e-6*6./(rho_bc*pi*d_bc**3)* &
                       exp(-9./2. * (log(sigma_bc))**2), 1.e10)
          ENDIF
          nbc = max(nbc, 1.e-10) 
          call S_max(Sat_max, nbc, w1, tc)

          if( tc < -20.0 ) then
            Ni_dust = MIN(MAX(imass(7)/dust_surf, dust_frac_min),   &
                     dust_frac_max)*EXP(-0.639 + 0.1296*(max(100.*(rhid - 1. ), Sat_max)))
          else
            Ni_dust = MIN(MAX(imass(7)/dust_surf, dust_frac_min),   &
                     dust_frac_max)*EXP(-0.639 + 0.1296*(100.*(rhid - 1. )))
          endif
! <---h1g, 2012-06-30

        endif 

!-------------------------------------------------------------------------
!    convert from 1/l to to  m-3
!-------------------------------------------------------------------------
        Ni_dust = 1.e3*Ni_dust

!-------------------------------------------------------------------------
!    if temp > -5C, no ice nucleation occurs.
!-------------------------------------------------------------------------
      else
        do_hom = .false.
      endif               

!-------------------------------------------------------------------------
!    sum the total available ice nulei. limit the number to be between 0 
!    and 1.0e8.
!-------------------------------------------------------------------------
      Ni = MIN( MAX(Ni_sulf + Ni_dust + Ni_bc, 0.), 1.e8) 

!-------------------------------------------------------------------------
!    define output variable hom to be 1 when homogeneous nucleation has 
!    occurred. note that this occur if do_hom is true for any of the 
!    segments of the pdf integration. 
!-------------------------------------------------------------------------
      if (do_hom) then
        hom = 1.
      endif 

!-------------------------------------------------------------------------


END SUBROUTINE ice_nucl_k



!##########################################################################

SUBROUTINE fast ( W1, naer, Ni, tc  )   

REAL, INTENT (IN)     :: W1, naer, tc
REAL, INTENT (INOUT ) :: Ni


      REAL, PARAMETER :: a1   = 0.0231
      REAL, PARAMETER :: a2_w = -1.639
      REAL, PARAMETER :: a2_c = -6.045
      REAL, PARAMETER :: b1   = -0.008
      REAL, PARAMETER :: b2_w = -0.042
      REAL, PARAMETER :: b2_c = -0.112
      REAL, PARAMETER :: c1  = 0.0739
      REAL, PARAMETER :: c2 = 1.237
 
!-------------------------------------------------------------------------
!--local variables-----

      REAL :: x1 , x2

      IF ( tc .GT. -64. ) THEN
        x1 =  a2_w + b2_w * tc +  c2 * LOG ( W1 )
        x2 = a1 + b1 * tc + c1 *  LOG ( W1 )
      ELSE
        x1 =  a2_c + b2_c * tc +  c2 * LOG ( W1 )
        x2 = a1 + b1 * tc + c1 *  LOG ( W1 )
      END IF
      Ni = MIN(EXP(x1)*naer**x2 , naer)



END SUBROUTINE fast


!#########################################################################

SUBROUTINE slow ( W1, naer, Ni ,tc )   

REAL, INTENT (IN)     :: W1, naer, tc
REAL, INTENT (INOUT ) :: Ni


      REAL, PARAMETER :: a1s   = -0.3949
      REAL, PARAMETER :: a2s   = 1.282
      REAL, PARAMETER :: b1s   = -0.0156
      REAL, PARAMETER :: b2s   = 0.0111
      REAL, PARAMETER :: b3s   = 0.0217
      REAL, PARAMETER :: c1s   = 0.120
      REAL, PARAMETER :: c2s   = 2.312

  
      REAL :: x1 , x2

      x1 = a2s + (b2s + b3s * LOG( W1 ))  * tc + c2s * LOG (W1) 
      x2 = a1s + b1s * tc + c1s * LOG(W1)
      Ni = MIN( EXP(x1)*naer**x2, naer)


END SUBROUTINE slow

!########################################################################

SUBROUTINE  bc_het(Ni, nbc, w1, tc)

REAL, INTENT (INOUT ) :: Ni
REAL, INTENT (IN )    :: NBC, W1, TC

      REAL, PARAMETER :: a11 = 0.0263
      REAL, PARAMETER :: b11 = -0.008
      REAL, PARAMETER :: a12 = -0.0185
      REAL, PARAMETER :: b12 = -0.0468
      REAL, PARAMETER :: a21 = 2.758
      REAL, PARAMETER :: b21 = -0.2667
      REAL, PARAMETER :: a22 = 1.3221
      REAL, PARAMETER :: b22 =  -1.4588

      REAL :: a, b 
 
      a = exp((a21*log(w1) + a22) + (a11*log(w1) + a12)*tc)
      b = (b21*log(w1) + b22) + (b11*log(w1) + b12)*tc 
      Ni = MIN(a*nbc**b, nbc) 


END SUBROUTINE  bc_het

!------------------------------------------------------------------------

!########################################################################

SUBROUTINE  S_max(Sat_max, nbc, w1, tc)
REAL, INTENT (INOUT ) :: Sat_max
REAL, INTENT (IN )    :: nbc, w1, tc
real                  ::  A,  B,   C,   &
                          a1, a2, a3,   &
                          b1, b2, b3,   &
                          c1, c2, c3

         a1 = -0.2035*nbc**(-0.8854)
         a2 =  0.2725*nbc**(-0.415)
         a3 = -0.0069

         b1 = -24.759*nbc**(-0.8831)
         b2 = 29.893*nbc**(-0.4067)
         b3 = -0.672

         c1 = -732.36*nbc**(-0.8712)
         c2 = 822.49*nbc**(-0.3951)
         c3 = 6.702

         A = a1*w1*w1 + a2 * w1 +a3
         B = b1*w1*w1 + b2 * w1 +b3
         C = c1*w1*w1 + c2 * w1 +c3

         Sat_max = A * tc*tc + B * tc + C

END SUBROUTINE  S_max
!------------------------------------------------------------------------


END MODULE ice_nucl_mod
