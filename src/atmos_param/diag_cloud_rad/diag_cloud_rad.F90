MODULE DIAG_CLOUD_RAD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!CLOUD RADIATIVE PROPERTIES
!
!       May-Oct  1998 -> Sep 2000
!       Contact persons: Tony Gordon, Bill Stern (for modified code)
!                        Steve Klein for original Fotran 90 code)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This module solves for the radiative properties of
!       every cloud.  In particular it uses either the
!       two stream approach or the delta-Eddington approach
!       to solve for the longwave emissivity, the ultra-violet-
!       visible reflectivities and absorptions, and the
!       near-infrared reflectivities and absorptions.
!
!       Modifications to Steve Klein's version have been made
!       to accomodate the empirical diagnostic cloud scheme of Tony Gordon,
!       frozen version v197 as discussed below.
!     
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------

use       mpp_mod, only: input_nml_file
use       fms_mod, only: file_exist, check_nml_error, open_namelist_file, &
                         mpp_pe, mpp_root_pe, close_file, &
                         write_version_number, stdlog

! Steve Klein's Cloud_Rad module
use Cloud_Rad_Mod, ONLY: CLOUD_RAD, CLOUD_RAD_INIT, cloud_rad_k_diag


!-------------------------------------------------------------------
 
  implicit none

!-------------------------------------------------------------------

!        The module contains the following:
!
!SUBROUTINES
!
! CLOUD_TAU_DRIVER
!   calls a sequence of suboutines to compute
!   cloud optical and radiative properties, as detailed
!   below -->
! CLOUD_PRES_THICK_FOR_TAU
!   computes cloud-type dependent set of pressure 
!   thicknesses for each distinct cloud layer, which are  
!   used to parameterize cloud optical depths
! CLOUD_OPTICAL_DEPTHS
!   Specify / crudely parameterize cloud optical depths
!   for distinct cloud layers,incorporating a 
!   parameterization scheme for non-anvil cirrus
!   proposed by Harshvardhan, based upon observations by
!   Platt and Harshvardhan
! CLOUD_OPTICAL_PROP_tg
!   for each cloud, it first establishes the standard
!   values of cloud optical depths (tau) in the 
!   vis+nir band for low, middle and high (anvil or non-
!   anvil cirrus) clouds. Then, it calculates the 
!   effective cloud liquid and cloud ice water paths, from
!   the cloud optical depths, mass aborption coefficients
!   for cloud water and cloud ice, and the temperature-
!   dependent ratio of cloud ice to cloud water path.
!   In the vis band, the single scattering albedo (wo), and
!   the asymmetry parameter (g) are both specified as 
!   constants (wo = 0.99999 and g = 0.85); 
!   In the nir band, wo is specified as a constant, 0.9942
!   (standard case), but is a function of zonal mean water 
!   vapor mixing ratio at the model's lowest vertical level
!   (anomalous absorption case). This subroutine also 
!   computes the longwave emissivity of each cloud.
! CLOUD_RAD   
!   solves for the radiative properties of the
!   every cloud.  In particular it uses either the
!   two stream approach or the delta-Eddington approach
!   to solve for the longwave emmissivity, the 
!   ultra-violet - visible reflectivities and absorptions,
!   and the near-infrared reflectivities and absorptions.
!   ****    Note    ****  This subroutine is from
!   the CLOUD_RAD_MOD module of Steve Klein .
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     PARAMETERS OF THE SCHEME
!
!     taumin       minimum permissible tau
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! 
!       PARAMETERS In Tony Gordon's v197 scheme only
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate).
!                       The diffusivity factor, 1.66 is incoporated into
!                       k_lw_liq.
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate).
!                       The diffusivity factor, 1.66 is incoporated into
!                       k_lw_ice.
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate).
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate).
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!       qsat_min      For zonal mean saturation water vapor mixing ratios
!                        less than this value entering into the variable 
!                        anomalous absorption calculations, the nir single 
!                        scattering albedo will asymptote to w0_anom2_nir
!                        where w0_anom2_nir is defined below.
!       qsat_trans    Transition value of zonal mean saturation mixing ratio
!                        between two branches of a piecewise continuous
!                        function entering into the variable anomalous
!                        absorption single scattering albedo calculations.
!       qsat_max      For zonal mean saturation water vapor mixing ratios
!                        greater than this value entering into the variable
!                        anomalous absorption calculations, the nir single
!                        scattering albedo will asymptote to w0_anom1_nir,
!                        where w0_anom1_nir is defined below.
!
!       w0_norm_uv    Normal, constant single scattering albedo in the uv-vis
!                        wavelength band of the radiation spectrum.
!       w0_norm_nir   Normal, constant single scattering albedo in the nir
!                        wavelength band of the radiation spectrum.
!       w0_anom1_nir  Asymptotic minimum value of single scattering albedo
!                        for variable anomalous absorption; also the
!                        low constant value, if L_anom_abs_g is set to TRUE.
!       w0_anom2_nir  Asymptotic maximum value of single scattering albedo
!                        for variable anomalous absorption,usually occurring
!                        at high latitudes.
! 
!       g_norm        Normal, constant asymmetry parameter used in Tony Gordon's
!                        v197 scheme. It is independent of wavelength.                                
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!          

 private



!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: diag_cloud_rad.F90,v 19.0 2012/01/06 20:05:09 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------

! REAL, PARAMETER :: taumin = 1.E-06
! Allow optical depths = 0.0 to remain 0.0
REAL, PARAMETER :: taumin = 0.0
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!     PARAMETERS in Tony Gordon's v197 scheme only:

!       Absorption Coefficients
REAL, PARAMETER     :: k_lw_liq = 140., k_lw_ice = 100.
REAL, PARAMETER     :: k_sw_liq = 130., k_sw_ice =  74.

!       Single scattering albedo and asymmetry parameters
REAL, PARAMETER     :: w0_norm_uv   = 0.99999, w0_norm_nir  = 0.9942
REAL, PARAMETER     :: w0_anom1_nir = 0.9700,  w0_anom2_nir = 0.9980
REAL, PARAMETER     :: g_norm = 0.85

!       Parameters controlling the proportion of ice to liquid phase
REAL, PARAMETER     :: tk_all_ice = 258.16, tk_all_liq = 268.16
 
!       Mixing ratio parameters in dimensionless units.
!       (To convert these parameters to g/kg, one would multiply them by 1000).

REAL, PARAMETER  :: qsat_min = 0.0, qsat_trans = 0.01, qsat_max = 0.02


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  An initialization subroutine, DIAG_CLOUD_RAD_INIT is called from 
!  clouds_tg_init do_crad_init will be reset to .false.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

!---------------------------------------------------------------------
! --- NAMELIST (clouds_rad_tg_nml)
!---------------------------------------------------------------------
!     l_har_anvil - logical variable = t -> anvil and super anvil cirrus
!                 clouds are treated as warm clouds.
!     l_har_coldcld - logical variable = t -> activates cold cloud portion
!                 of Harshvardhan scheme.
!                 Default value has been changed to .true., consistent with
!                 AMIP 2 run and preferred when employing RAS with no cap.
!                 l_har_coldcld = .false. would correspond to
!                 v197 frozen model runs.
!     l_anom_abs_g - logical variable = t -> anomalous absorption is
!                  specified as a constant value of single scattering albedo.
!     l_anom_abs_v - logical variable = t -> anomalous absorption is
!                  computed as a piecewise continuous function of zonal
!                  mean saturation mixing ratio.
!---------------------------------------------------------------------

 logical :: &
     & l_har_anvil = .true.,l_har_coldcld = .true., &
     & l_anom_abs_g = .false., l_anom_abs_v = .false.

  NAMELIST / diag_cloud_rad_nml /  &                          
     & l_har_anvil, l_har_coldcld, l_anom_abs_g, l_anom_abs_v

!  Need to create an initialization subroutine which reads this namelist.
!  Also, in that subroutine, check that L_anom_abs_g and L_anom_abs_v
!  are not both set to .true., after reading the namelist.
!  If they are, either code an ERROR EXIT, or else, allow L_anom_abs_v to take
!  precedence, i.e., reset L_anom_abs_g = .false. :

!     IF (L_anom_abs_v) THEN
!        L_anom_abs_g = .FALSE.
!    ENDIF   

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!    ******    For Tony Gordon's v197 scheme only:    ****
!
!    The maximum number of distinct cloud layers in a vertical column,
!    max_cld = maxval(nclds), can be supplied by the cldtim driver in module
!    clouds_tg. Also, the cldtim driver can call a special subroutine to
!    compute the useful diagnostic total_cld_amt.
!
!    For efficiency, it is preferrable to work with the compressed vertical
!    index k'. Also, all of the argument calls in this module are presently
!    set up that way. Of course, that can be changed, by filling all levels
!    within a distinct cloud layer with the same values of data
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 public cloud_tau_driver,diag_cloud_rad_init, diag_cloud_rad_end, &
     cloud_pres_thick_for_tau, cloud_optical_depths,  &
     cloud_optical_depths2, &
     cloud_opt_prop_tg_lw, cloud_opt_prop_tg_sw, cloud_opt_prop_tg2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 contains


!########################################################################
!########################################################################

SUBROUTINE CLOUD_TAU_DRIVER (qmix_kx, tempcld, tau, coszen,  &
                             r_uv, r_nir, ab_uv, ab_nir, em_lw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following radiative properties of clouds
!
!               1. r_uv:   cloud reflectance in uv band
!               2. r_nir:  cloud reflectance in nir band
!               3. ab_uv:  cloud absorption in uv band
!               4. ab_nir: cloud absorption in nir band
!               5. em_lw:  longwave cloud emissivity
!
!               Note: Our nir is split, later, into 3 smaller bands
!                     employed in a module created by Steve Klein.
!                     The radiative properties do not vary amongst
!                     those bands.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!      NCLDS       number of (random overlapping) clouds in column
!      ICLD          marker array of cloud types/heights (at cloud levels)
!      CLDTOP     indices of model levels where cloud tops of distinct cloud
!                    layers are located. 
!      CLDBAS     indices of model levels where cloud bases of distinct cloud
!                    layers are located.
!      CLDAMT       cloud amount fraction of distinct cloud layers
!      DELP_TRUE  true cloud pressure thickness of distinct cloud layers
!      TEMPCLD    cloud layer mean temperature (degrees Kelvin) of distinct
!                    cloud layers
!      PFULL        pressure at full levels (Pascals)
!      PHALF        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!      PSFC       Surface pressure field
!      COSZEN     cosine of the zenith angle

! -----------------------------------------------------------------------------
! as of May 1999 these variables are internal to Steve Klein's routine

!      l2strem      logical variable indicating 2 stream operating or not
!                          l2strem = T  2 stream solution to the calculation
!                          of cloud raditive properties
!                          l2strem = F  Delta-Eddington solution to the
!                          calculation of cloud radiative properties.
!
!                            IF l2strem = T then the solution does not
!                            depend on solar zenith angle
!
!                    [ namelist variable in Steve Klein's Cloud_Rad module ]
!
!      taucrit      critical tau for switching direct beam to diffuse beam
!
!                    [ namelist variable in Steve Klein's Cloud_Rad module ]
! -----------------------------------------------------------------------------
!
!       -------------
!OUTPUT:
!       -------------
!
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!       em_lw        longwave cloud emissivity
!
!       -------------
!       The following variables might be elevated to "OUTPUT" status later.
!          They are computed in subroutines called by this subroutine.
!          Perhaps they would be needed for diagnostics purposes.
!          If so, add w0 and gg to the argument list of CLOUD_TAU_DRIVER
!
!       -------------
!       w0           single scattering albedo for each band
!       gg           asymmetry parameter for each band
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t        looping variable
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of vertical levels
!       max_cld      maximum number of distinct cloud layers in whole array
!       LWP          cloud liquid water path (kg of condensate per square meter)
!       IWP          cloud ice water path (kg of condensate per squre meter)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!      For Tony Gordon's v197 scheme only
!
!      L_anom_abs-g Logical namelist variable. If true, anomalous absorption
!                   is represented by a constant value of single scattering
!                   albedo. The default value and the v197 setting are
!                   both false.
!
!      L_anom_abs_v Logical namelist variable. If true, anomalous absorption
!                   is computed as a piecewise continuous function of zonal
!                   mean saturation water vapor mixing ratio at the model's
!                   vertical level closest to the earth's surface. The
!                   default value is false. The analogous namelist variable
!                   in v197, LWALBV is set to TRUE. 
!                   If both L_anom_abs_g and L_anom_abs_v are set to FALSE,
!                   then the single scattering albedo, w0,  assumes its normal
!                   constant value of 0.9942, in tg's version. w0 is a
!                   variable in Steve Klein's version.
!                   L_anom_abs_v takes precedence over L_anom_abs_g, if
!                   the former is TRUE and the latter is FALSE.
!
!     qmix_kx       water vapor mixing ratio at the
!                   model's vertical level closest to the earth's surface.
!                   It is more convenient to pass the zonal mean, in case
!                   single column tests are performed. Alternatively,
!                   the zonal mean could be computed within this subroutine. 
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

real, intent(in), dimension(:,:,:) ::    tempcld
real,    intent (in),dimension(:,:)   :: qmix_kx

real,    intent(in), dimension (:,:)  ::  coszen


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real,     intent (out),dimension(:,:,:), optional  :: r_uv,r_nir,  &
                           ab_uv,ab_nir,em_lw

!  *****************************************************************

real,     intent (in   ),dimension(:,:,:,:) :: tau

real, dimension(size(tau,1),size(tau,2),size(tau,3),size(tau,4)) :: w0,gg

real, dimension(size(tau,1),size(tau,2),size(tau,3)) :: lwp,iwp

!  *****************************************************************




!compute cloud radiative properties



if (present(em_lw)) then
          call cloud_opt_prop_tg(lwp,iwp,tau,w0,gg,tempcld,qmix_kx,  &
       em_lw=em_lw)
        else
          call cloud_opt_prop_tg(lwp,iwp,tau,w0,gg,tempcld,qmix_kx)
endif
 
!  From Steve Klein's cloud_rad_mod

        if (present(r_uv)) then
            call cloud_rad(tau,w0,gg,coszen,r_uv,r_nir,ab_uv,ab_nir)
        endif

!--------------------------------------------------------------------


END SUBROUTINE CLOUD_TAU_DRIVER





!########################################################################
subroutine cloud_pres_thick_for_tau (nclds,icld,cldtop,cldbas, &
     &          delp_true,lhight,lhighb, lmidt, lmidb, llowt,lk,delp, &
     &          phalf, psfc )

! This subroutine calculates a special cloud-type dependent set of 
! cloud pressure thicknesses that will be employed to compute cloud optical
! depths.

!===================================================================

!  parameters used in this routine 
!
! standard optical pressure thickness quantities (hPa)
REAL, PARAMETER :: delp_high = 31.25, delp_mid = 75.0, delp_low = 112.5
REAL, PARAMETER :: delp_high_thk = 2.0*delp_high, delp_cnv = 112.5
REAL, PARAMETER :: delp_high_thk_anvil = 2.0*delp_high
! REAL, PARAMETER :: delp_strat = 112.5, delp_min = 25.0, delp_thk_crit = 125.0
! impose v197 value for delp_min
REAL, PARAMETER :: delp_strat = 112.5, delp_min = 9.0, delp_thk_crit = 125.0

! delp_high - standard specified pressure thickness for high (warm) clouds,
!             which is used to compute their optical depth.
! delp_mid -  standard specified pressure thickness for middle (warm) clouds,
!             which is used to compute their optical depth.
! delp_low -  standard specified pressure thickness for low (warm) clouds, 
!             which is used to compute their optical depth, 
!             provided that those clouds are not too close to the ground.
! delp_high_thk - specified pressure thickness  for high clouds (other than 
!             precipitating convective clouds), whose true cloud pressure
!             thickness exceeds delp_thk_crit, and 
!             which is used to compute their optical depth.
! delp_high_thk_anvil - specified pressure thickness  for anvil cirrus clouds 
!             whose true cloud pressure thickness exceeds delp_thk_crit, and
!             which is used to compute their optical depth.
! delp_cnv -  cloud optical depth for precipitating convective clouds with
!             bases in the "low" cloud region of the atmosphere,
!             which is used to compute their optical depth.
! delp_strat - standard specified pressure thickness for marine stratus clouds,
!             which is used to compute their optical depth.
!             They could be quite optically thick, in some regions, e.g.,
!             off the coast of Peru, based upon ISCCP obs.
! delp_min -  minimum pressure thickness of a cloud layer.
! delp_thk_crit - critical true cloud pressure thickness, which, if exceeded,
!             signifies a thick stratiform cloud.

!===================================================================


! Arguments (intent in)

real,    intent(in), dimension(:,:,:) :: phalf, delp_true
integer, intent(in), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(in), dimension(:,:)   ::  nclds
integer, intent(in), dimension(:,:)   :: lhight, lhighb, lmidt, lmidb, llowt
integer, intent(in)                   :: lk
real,    intent(in), dimension (:,:)  ::  psfc

!      INPUT
!      ------

!       PHALF        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       PSFC         Surface pressure field
!       NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!       ICLD          marker array of cloud types/heights (at cloud levels)
!       CLDTOP     index of cloud tops (at cloud levels)
!       CLDBAS     index of cloud bottoms (at cloud levels)
!       DELP_TRUE  true cloud pressure thickness of distinct cloud layers 
!                    (at cloud levels - in Pa)
!       LHIGHT        vertical level index upper limit for high cloud tops
!       LHIGHB        vertical level index lower limit for high cloud bases
!       LMIDT         vertical level index upper limit for mid cloud tops
!       LMIDB         vertical level index lower limit for mid cloud bases
!       LLOWT         vertical level index upper limit for low cloud tops
!       LK            vertical level below which no low cloud bases can exist

!===================================================================

! Arguments (intent out)

real, intent(out), dimension(:,:,:) :: delp

!      OUTPUT
!      ------

!      DELP     cloud pressure thickness used to calculate cloud optical depths 
!                   of distinct cloud layers
!=======================================================================
!  (Intent local)

 integer, dimension (size(cldtop,1),size(cldtop,2)) :: cldt,cldb,icwork
 real, dimension (size(cldtop,1),size(cldtop,2)) :: delpstd,pwork

! scalars
 integer k, max_cld, idim, jdim, kmaxp

 integer i, j

!      DELP_STD     tentative cloud pressure thickness values 
!                   (at cloud levels)
!      CLDT, CLDB   work arrays for cloud top,base indices
!      ICWORK       cloud type marker work array at cloud levels 
!      PWORK        work array for pressure values at cloud levels (hPa)

!===================================================================

! define horizontal dimensions

  idim = SIZE( cldtop, 1 )
  jdim = SIZE( cldtop, 2 )
  kmaxp = SIZE( phalf, 3 )


! Initialize cloud pressure thickness array
      delp = 0.0

! find maximum number of cloud levels
      max_cld  = maxval(nclds(:,:))
                   if (max_cld .ge. 1) then

!-----------------------------------------------------------------------
! <><><><><><><><>   calc cloud press thickness <><><><><><><><>
!-----------------------------------------------------------------------

      do k = 1,max_cld

! Initialize internal arrays
      cldt = kmaxp
      cldb = 0
      icwork = 0
      pwork = 0.0
      delpstd = 0.0
         where (nclds(:,:) .ge. k)
           cldt(:,:) = cldtop(:,:,k)      
           cldb(:,:) = cldbas(:,:,k) 
         end where

! separate cloud type marker from cloud height marker
         icwork(:,:) = mod(icld(:,:,k),100)

! fill array pwork with the pressure difference value of the half level 
! corresponding to the cloud top level cldt and the surface pressure
! also convert from Pa to hPa
         do j=1,jdim
         do i=1,idim
           pwork(i,j) = (psfc(i,j) - phalf(i,j,cldt(i,j)))*.01
         end do
         end do

         where (cldb(:,:) .ge. lhight(:,:) .and. cldb(:,:) .le. lhighb(:,:) )
           delpstd(:,:) = delp_high
         end where

         where (cldb(:,:) .ge. lmidt(:,:) .and. cldb(:,:) .le. lmidb(:,:) )
           delpstd(:,:) = delp_mid
         end where

         where ( (cldb(:,:) .ge. llowt(:,:) .and. cldb(:,:) .le. lk) .and. &
     &                 (pwork(:,:) .ge. delp_low) )
           delpstd(:,:) = delp_low
         end where

         where ( (cldb(:,:) .ge. llowt(:,:) .and. cldb(:,:) .le. lk) .and. &
     &                 (pwork(:,:) .lt. delp_low) )
           delpstd(:,:) = max(pwork(:,:), delp_min)
         end where


! tentative value of cloud pressure thickness used in cloud optical depth calc
         where (cldb(:,:) .ge. lhight(:,:) .and. cldb(:,:) .le. lk)
           delp(:,:,k) = delpstd(:,:)
         end where
! redefine cloud pressure thickness for precipitating convective clouds with
! low cloud bases so that a small value such as delp_min is not a possiblity
         where (cldb(:,:) .ge. llowt(:,:) .and. cldb(:,:) .le. (lk-1) .and. &
     &               icwork(:,:) .eq. 5)
           delp(:,:,k) = delp_cnv
         end where

! use pwork to store delp_true  converted from Pa to hPa
         pwork(:,:) = delp_true(:,:,k)*.01

! redefine cloud pressure thickness for thick, high stratiform clouds
         where ( (cldb(:,:) .ge. lhight(:,:) .and. &
     &            cldb(:,:) .le. lhighb(:,:)) .and. &
     &           (pwork(:,:) .gt. delp_thk_crit) .and. &
     &           (icwork(:,:) .eq. 1) )
           delp(:,:,k) = delp_high_thk
         end where
! redefine cloud pressure thickness for thick, anvil cirrus clouds
         where ( (cldb(:,:) .ge. lhight(:,:) .and. &
     &            cldb(:,:) .le. lhighb(:,:)) .and. &
     &           (pwork(:,:) .gt. delp_thk_crit) .and. &
     &           (icwork(:,:) .ge. 6) )
           delp(:,:,k) = delp_high_thk_anvil
         end where
! redefine cloud pressure thickness for marine stratus clouds
         where ( (cldb(:,:) .ge. llowt(:,:) .and. cldb(:,:) .le. lk) .and. &
     &           (cldb(:,:) .eq. cldt(:,:)) .and. (icwork(:,:) .eq. 3) )
           delp(:,:,k) = min(delpstd(:,:),delp_strat)
         end where

      end do

      
                   endif

end subroutine cloud_pres_thick_for_tau

!########################################################################

subroutine cloud_optical_depths2 (nclds,icld,cldtop,cldbas,tempcld,delp, &
     &          tau,phalf,           liq_frac )

! This subroutine specifies/crudely parameterizes cloud optical depths 
! of non-anvil cirrus clouds based upon a parameterization scheme 
! resembling that proposed by Harshvardhan,
! (which is based on observations by Platt and Harshvardhan).

!===================================================================

!  namelist quantities used in this routine (defined in module intro above).

!     l_har_anvil - logical variable = t -> anvil and super anvil cirrus
!                 clouds are treated as warm clouds.
!     l_har_coldcld - logical variable = t -> activates cold cloud portion
!                 of Harshvardhan scheme.
!===================================================================

!  parameters used in this routine 

REAL, PARAMETER :: ctok = 273.16, t_ref = ctok - 82.5 
REAL, PARAMETER :: t_warm = ctok, t_cold = ctok-10.0
REAL, PARAMETER :: temp_dif_min = 1.0, inv_delta_t = 1./(t_warm-t_cold)
REAL, PARAMETER :: harshb_std = 0.08, harshb_cnv = 2.0*harshb_std
REAL, PARAMETER :: harshb_anvil = harshb_std, harshb_super_anvil = harshb_std
! The following value is 6 times larger than suggested by Platt-Harshvardan.
! The new value boosts the cold cloud emissivities in the tropics, giving
! closer agreement with ERBE.
REAL, PARAMETER :: harsha_cold = 12.8e-06

! ctok -      freezing point in deg. K..
! t_ref -     a reference minimum permitted cloud temperature in the
!             Harshvardan parameterization.
! t_warm -    the lowest temperature the cloud layer can have and still 
!             be considered warm.
! t_cold -    the highest temperature the cloud layer can have and still 
!             be considered cold. 
! temp_dif_min - the minimum temperature difference to be used in the 
!             computation of Platt-Harshvardhan cold cloud optical depths.
! inv_delta_t =  1./(t_warm-t_cold)
! harshb_std - Harshvardan coefficient for a "standard" cloud layer 
! harshb_cnv - Harshvardan coefficient for convective clouds treated as warm 
! harshb_anvil - Harshvardan coefficient for anvil cirrus clouds 
! harshb_super_anvil - Harshvardan coefficient for super anvil cirrus clouds 
! harsha_cold - Harshvardan coefficient for a cold high cloud 

!===================================================================


! Arguments (intent in)

real, intent(in), dimension(:,:,:)    :: phalf, delp, tempcld
integer, intent(in), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(in), dimension(:,:)   ::  nclds

!      INPUT
!      ------

!       PHALF        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!      ICLD          marker array of cloud types (at cloud levels)
!      CLDTOP     index of cloud tops (at cloud levels)
!      CLDBAS     index of cloud bottoms (at cloud levels)
!      TEMPCLD    cloud layer mean temperature 
!                 (degrees Kelvin, at cloud levels)
!      DELP      cloud pressure thickness used for cloud optical depth 
!                   (at cloud levels)

!===================================================================

! Arguments (intent out)

real, intent(inout), dimension(:,:,:,:) :: tau
real, intent(inout), dimension(:,:,:) ::           liq_frac

!      OUTPUT
!      ------

!      TAU     cloud optical depth (at cloud levels)

!=======================================================================

!  (Intent local)

 integer, dimension (size(cldtop,1),size(cldtop,2)) :: cldt,cldb,icwork
 real, dimension (size(cldtop,1),size(cldtop,2)) :: & 
                                          temp_dif,tau_cold,tau_warm,pwork
 real, dimension (size(tau,1),size(tau,2),size(tau,3) )  :: tau_vis

! scalars
 integer i, j, idim, jdim, kmaxp
 integer k, max_cld
 integer n, max_band

! real 

!      CLDT, CLDB   work arrays for cloud top,base indices
!      TAU_COLD     cloud optical depth work array for cold calculation in
!                   transition regime 
!      TAU_WARM     cloud optical depth work array for warm calculation in
!                   transition regime 
!      TEMP_DIF     work array for temperature diference of distinct cloud
!                   layers 
!      ICWORK       cloud type marker work array of distinct cloud layers
!      PWORK        work array for pressure values at cloud levels (hPa) 
!      TAU_VIS      work array for cloud optical depth in visible part
!                   of spectrum
!      MAX_BAND     maximum number of radiative bands for cloud optical depth
!      MAX_CLD      maximum number of distinct cloud layers within a 
!                   vertical column

!===================================================================


! define horizontal dimensions

  idim = SIZE( cldtop, 1 )
  jdim = SIZE( cldtop, 2 )
  kmaxp = SIZE( phalf, 3 )


! find maximum number of cloud levels
      max_cld  = maxval(nclds(:,:))

! define maximum number of wave number bands for tau
      max_band = size(tau,4)


!===================================================================

! Initialize working and final cloud optical depth arrays
      tau_vis = 0.0
      tau     = 0.0

! find maximum number of clouds
      max_cld  = maxval(nclds(:,:))
                   if (max_cld .ge. 1) then

!-----------------------------------------------------------------------
! <><><><><><><><>   calc cloud optical depths <><><><><><><><>
!-----------------------------------------------------------------------

      do k = 1,max_cld

! Initialize internal arrays
      cldt = kmaxp
      cldb = 0
      icwork = 0
         where (nclds(:,:) .ge. k)
           cldt(:,:) = cldtop(:,:,k)      
           cldb(:,:) = cldbas(:,:,k) 
         end where

! fill array pwork with pressure values at the half level corresponding
! to the cloud top level cldt and convert from Pa to hPa
         do j=1,jdim
         do i=1,idim
           pwork(i,j) = phalf(i,j,cldt(i,j))*.01
         end do
         end do

! separate cloud type marker from cloud height marker
         icwork(:,:) = mod(icld(:,:,k),100)

! Standard case: warm cloud treatment will be applied to warm and cold clouds.
! Later, the cold cloud cases will be re-computed, if namelist parameter
! l_har_coldcld is set to true.

! preliminary optical depths computed everywhere

         tau_vis(:,:,k) = harshb_std * delp(:,:,k)

! redefine tau_vis for convective clouds, which are always treated as warm
         where (icwork(:,:) .eq. 5)
           tau_vis(:,:,k) = harshb_anvil * delp(:,:,k)
         end where

                   if (l_har_anvil) then
! redefine tau_vis for high clouds that meet the anvil cirrus criterion
         where (icld(:,:,k) .eq. 106)
           tau_vis(:,:,k) = harshb_anvil * delp(:,:,k)
         end where

! redefine tau_vis for high clouds that meet the super anvil cirrus criterion
         where (icld(:,:,k) .eq. 107)
           tau_vis(:,:,k) = harshb_super_anvil * delp(:,:,k)
         end where
                   endif

                   if (l_har_coldcld .and. l_har_anvil ) then

! ordinary rhum cirrus clouds: redefine tau_vis for cold regime and
! transition regime
         where (icld(:,:,k) .eq. 101 .and. (tempcld(:,:,k).le.t_cold))
! compute Platt-Harshvardhan optical depths in cold cloud regime
           temp_dif(:,:) = max ( (tempcld(:,:,k)-t_ref), temp_dif_min)
           tau_vis(:,:,k) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
         end where

         where (icld(:,:,k) .eq. 101 .and.  &
     &             (tempcld(:,:,k).gt.t_cold .and. tempcld(:,:,k).lt.t_warm))
! compute Platt-Harshvardhan optical depths in the transition regime,
! by linearly interpolating solutions from the cold and warm regimes
! with respect to temperature
           temp_dif(:,:) = t_cold - t_ref
           tau_cold(:,:) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
           tau_warm(:,:) = harshb_std * delp(:,:,k)
           tau_vis(:,:,k) = tau_cold(:,:) + (tau_warm(:,:)-tau_cold(:,:)) * &
     &              (tempcld(:,:,k)-t_cold) * inv_delta_t 
         end where

                   endif

                   if (l_har_coldcld .and. (.not.l_har_anvil) ) then

! ordinary rhum cirrus clouds: redefine tau_vis for cold regime and
! transition regime
         where ( (icld(:,:,k) .eq. 101 .or. icld(:,:,k) .eq. 106 .or. &
     &            icld(:,:,k) .eq. 107) .and. &
     &              (tempcld(:,:,k).le.t_cold))
! compute Platt-Harshvardhan optical depths in cold cloud regime
           temp_dif(:,:) = max ( (tempcld(:,:,k)-t_ref), temp_dif_min)
           tau_vis(:,:,k) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
         end where

         where ( (icld(:,:,k) .eq. 101 .or. icld(:,:,k) .eq. 106 .or. &
     &            icld(:,:,k) .eq. 107) .and. &
     &             (tempcld(:,:,k).gt.t_cold .and. tempcld(:,:,k).lt.t_warm))
! compute Platt-Harshvardhan optical depths in the transition regime,
! by linearly interpolating solutions from the cold and warm regimes
! with respect to temperature
           temp_dif(:,:) = t_cold - t_ref
           tau_cold(:,:) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
           tau_warm(:,:) = harshb_std * delp(:,:,k)
           tau_vis(:,:,k) = tau_cold(:,:) + (tau_warm(:,:)-tau_cold(:,:)) * &
     &              (tempcld(:,:,k)-t_cold) * inv_delta_t 
         end where

                   endif


      end do


!  The cloud optical depths of all four wave number bands are set equal to 
!  tau_vis, i.e., the value in the visible band.

     do n=1,max_band
           tau(:,:,:,n) = tau_vis(:,:,:)
     end do

        WHERE (tau(:,:,:,:) .lt. taumin)
               tau(:,:,:,:) = taumin
        END WHERE

         call  CLOUD_OPT_PROP_tg2 (tau, tempcld,           liq_frac)

                   else
     liq_frac = 0.0
                   endif

end subroutine cloud_optical_depths2




subroutine cloud_optical_depths (nclds,icld,cldtop,cldbas,tempcld,delp, &
     &          tau,phalf )

! This subroutine specifies/crudely parameterizes cloud optical depths 
! of non-anvil cirrus clouds based upon a parameterization scheme 
! resembling that proposed by Harshvardhan,
! (which is based on observations by Platt and Harshvardhan).

!===================================================================

!  namelist quantities used in this routine (defined in module intro above).

!     l_har_anvil - logical variable = t -> anvil and super anvil cirrus
!                 clouds are treated as warm clouds.
!     l_har_coldcld - logical variable = t -> activates cold cloud portion
!                 of Harshvardhan scheme.
!===================================================================

!  parameters used in this routine 

REAL, PARAMETER :: ctok = 273.16, t_ref = ctok - 82.5 
REAL, PARAMETER :: t_warm = ctok, t_cold = ctok-10.0
REAL, PARAMETER :: temp_dif_min = 1.0, inv_delta_t = 1./(t_warm-t_cold)
REAL, PARAMETER :: harshb_std = 0.08, harshb_cnv = 2.0*harshb_std
REAL, PARAMETER :: harshb_anvil = harshb_std, harshb_super_anvil = harshb_std
! The following value is 6 times larger than suggested by Platt-Harshvardan.
! The new value boosts the cold cloud emissivities in the tropics, giving
! closer agreement with ERBE.
REAL, PARAMETER :: harsha_cold = 12.8e-06

! ctok -      freezing point in deg. K..
! t_ref -     a reference minimum permitted cloud temperature in the
!             Harshvardan parameterization.
! t_warm -    the lowest temperature the cloud layer can have and still 
!             be considered warm.
! t_cold -    the highest temperature the cloud layer can have and still 
!             be considered cold. 
! temp_dif_min - the minimum temperature difference to be used in the 
!             computation of Platt-Harshvardhan cold cloud optical depths.
! inv_delta_t =  1./(t_warm-t_cold)
! harshb_std - Harshvardan coefficient for a "standard" cloud layer 
! harshb_cnv - Harshvardan coefficient for convective clouds treated as warm 
! harshb_anvil - Harshvardan coefficient for anvil cirrus clouds 
! harshb_super_anvil - Harshvardan coefficient for super anvil cirrus clouds 
! harsha_cold - Harshvardan coefficient for a cold high cloud 

!===================================================================


! Arguments (intent in)

real, intent(in), dimension(:,:,:)    :: phalf, delp, tempcld
integer, intent(in), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(in), dimension(:,:)   ::  nclds

!      INPUT
!      ------

!       PHALF        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!      ICLD          marker array of cloud types (at cloud levels)
!      CLDTOP     index of cloud tops (at cloud levels)
!      CLDBAS     index of cloud bottoms (at cloud levels)
!      TEMPCLD    cloud layer mean temperature 
!                 (degrees Kelvin, at cloud levels)
!      DELP      cloud pressure thickness used for cloud optical depth 
!                   (at cloud levels)

!===================================================================

! Arguments (intent out)

real, intent(out), dimension(:,:,:,:) :: tau

!      OUTPUT
!      ------

!      TAU     cloud optical depth (at cloud levels)

!=======================================================================

!  (Intent local)

 integer, dimension (size(cldtop,1),size(cldtop,2)) :: cldt,cldb,icwork
 real, dimension (size(cldtop,1),size(cldtop,2)) :: & 
                                          temp_dif,tau_cold,tau_warm,pwork
 real, dimension (size(tau,1),size(tau,2),size(tau,3) )  :: tau_vis

! scalars
 integer i, j, idim, jdim, kmaxp
 integer k, max_cld
 integer n, max_band

! real 

!      CLDT, CLDB   work arrays for cloud top,base indices
!      TAU_COLD     cloud optical depth work array for cold calculation in
!                   transition regime 
!      TAU_WARM     cloud optical depth work array for warm calculation in
!                   transition regime 
!      TEMP_DIF     work array for temperature diference of distinct cloud
!                   layers 
!      ICWORK       cloud type marker work array of distinct cloud layers
!      PWORK        work array for pressure values at cloud levels (hPa) 
!      TAU_VIS      work array for cloud optical depth in visible part
!                   of spectrum
!      MAX_BAND     maximum number of radiative bands for cloud optical depth
!      MAX_CLD      maximum number of distinct cloud layers within a 
!                   vertical column

!===================================================================


! define horizontal dimensions

  idim = SIZE( cldtop, 1 )
  jdim = SIZE( cldtop, 2 )
  kmaxp = SIZE( phalf, 3 )


! find maximum number of cloud levels
      max_cld  = maxval(nclds(:,:))

! define maximum number of wave number bands for tau
      max_band = size(tau,4)


!===================================================================

! Initialize working and final cloud optical depth arrays
      tau_vis = 0.0
      tau     = 0.0

! find maximum number of clouds
      max_cld  = maxval(nclds(:,:))
                   if (max_cld .ge. 1) then

!-----------------------------------------------------------------------
! <><><><><><><><>   calc cloud optical depths <><><><><><><><>
!-----------------------------------------------------------------------

      do k = 1,max_cld

! Initialize internal arrays
      cldt = kmaxp
      cldb = 0
      icwork = 0
         where (nclds(:,:) .ge. k)
           cldt(:,:) = cldtop(:,:,k)      
           cldb(:,:) = cldbas(:,:,k) 
         end where

! fill array pwork with pressure values at the half level corresponding
! to the cloud top level cldt and convert from Pa to hPa
         do j=1,jdim
         do i=1,idim
           pwork(i,j) = phalf(i,j,cldt(i,j))*.01
         end do
         end do

! separate cloud type marker from cloud height marker
         icwork(:,:) = mod(icld(:,:,k),100)

! Standard case: warm cloud treatment will be applied to warm and cold clouds.
! Later, the cold cloud cases will be re-computed, if namelist parameter
! l_har_coldcld is set to true.

! preliminary optical depths computed everywhere

         tau_vis(:,:,k) = harshb_std * delp(:,:,k)

! redefine tau_vis for convective clouds, which are always treated as warm
         where (icwork(:,:) .eq. 5)
           tau_vis(:,:,k) = harshb_anvil * delp(:,:,k)
         end where

                   if (l_har_anvil) then
! redefine tau_vis for high clouds that meet the anvil cirrus criterion
         where (icld(:,:,k) .eq. 106)
           tau_vis(:,:,k) = harshb_anvil * delp(:,:,k)
         end where

! redefine tau_vis for high clouds that meet the super anvil cirrus criterion
         where (icld(:,:,k) .eq. 107)
           tau_vis(:,:,k) = harshb_super_anvil * delp(:,:,k)
         end where
                   endif

                   if (l_har_coldcld .and. l_har_anvil ) then

! ordinary rhum cirrus clouds: redefine tau_vis for cold regime and
! transition regime
         where (icld(:,:,k) .eq. 101 .and. (tempcld(:,:,k).le.t_cold))
! compute Platt-Harshvardhan optical depths in cold cloud regime
           temp_dif(:,:) = max ( (tempcld(:,:,k)-t_ref), temp_dif_min)
           tau_vis(:,:,k) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
         end where

         where (icld(:,:,k) .eq. 101 .and.  &
     &             (tempcld(:,:,k).gt.t_cold .and. tempcld(:,:,k).lt.t_warm))
! compute Platt-Harshvardhan optical depths in the transition regime,
! by linearly interpolating solutions from the cold and warm regimes
! with respect to temperature
           temp_dif(:,:) = t_cold - t_ref
           tau_cold(:,:) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
           tau_warm(:,:) = harshb_std * delp(:,:,k)
           tau_vis(:,:,k) = tau_cold(:,:) + (tau_warm(:,:)-tau_cold(:,:)) * &
     &              (tempcld(:,:,k)-t_cold) * inv_delta_t 
         end where

                   endif

                   if (l_har_coldcld .and. (.not.l_har_anvil) ) then

! ordinary rhum cirrus clouds: redefine tau_vis for cold regime and
! transition regime
         where ( (icld(:,:,k) .eq. 101 .or. icld(:,:,k) .eq. 106 .or. &
     &            icld(:,:,k) .eq. 107) .and. &
     &              (tempcld(:,:,k).le.t_cold))
! compute Platt-Harshvardhan optical depths in cold cloud regime
           temp_dif(:,:) = max ( (tempcld(:,:,k)-t_ref), temp_dif_min)
           tau_vis(:,:,k) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
         end where

         where ( (icld(:,:,k) .eq. 101 .or. icld(:,:,k) .eq. 106 .or. &
     &            icld(:,:,k) .eq. 107) .and. &
     &             (tempcld(:,:,k).gt.t_cold .and. tempcld(:,:,k).lt.t_warm))
! compute Platt-Harshvardhan optical depths in the transition regime,
! by linearly interpolating solutions from the cold and warm regimes
! with respect to temperature
           temp_dif(:,:) = t_cold - t_ref
           tau_cold(:,:) = harsha_cold * temp_dif(:,:)**2 * delp(:,:,k)
           tau_warm(:,:) = harshb_std * delp(:,:,k)
           tau_vis(:,:,k) = tau_cold(:,:) + (tau_warm(:,:)-tau_cold(:,:)) * &
     &              (tempcld(:,:,k)-t_cold) * inv_delta_t 
         end where

                   endif


      end do


!  The cloud optical depths of all four wave number bands are set equal to 
!  tau_vis, i.e., the value in the visible band.

     do n=1,max_band
           tau(:,:,:,n) = tau_vis(:,:,:)
     end do

        WHERE (tau(:,:,:,:) .lt. taumin)
               tau(:,:,:,:) = taumin
        END WHERE
                   endif


end subroutine cloud_optical_depths

!########################################################################

subroutine CLOUD_OPT_PROP_tg(LWP,IWP,       &
                        tau,w0,gg,          &
                        tempcld,qmix_kx, em_lw )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau    :optical depth in each band
!               2. w0     :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw  :longwave cloud emissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
! *************************    WARNING    *****************************
!
!   The above bands are used by Steve Klein.
!   We retain the scheme from the v197 frozen model,instead.
!   Nominally, our band 2 is expanded into bands 2 + 3 + 4 of Slingo.
!   The same cloud optical depth is specified in all 4 bands.
!   The same asymmetry parameter is specified in all 4 bands.
!   For single scattering albedo w0, the uv value is specified in band 1,
!   while the nir value is specified in bands 2, 3, and 4.

!   ****  WARNING    ****  The code is intended to be applied to 2 to 4 bands.
!   An error check to check that this condition is satisfied is advised.  

! *********************************************************************
!            BAND               v197     
!             1               0.25-0.70
!             2               0.70-4.00
! *********************************************************************
!
!

!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = 0.9942  (in v197 - standard)
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = F(Z.M. mixing ratio at lowest model level)
!                                    (in v197 - anomalous absorption)
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!           g(:,:,1:2) = 0.85 in v197
!   
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation could be rewritten, after algebraic manipulation, as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    However, the other form of the equation, i.e., 
!    1 - exp(tau_liq + tau_ice) will actually be solved.

! *******************************************************************
!
!
!   (6)  v197 only: Must first solve for LWP and IWP knowing
!                   tau, k_sw_liq, k_sw_ice, wgt_liq and wgt_ice.
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!


!      L_anom_abs-g Logical namelist variable. If true, anomalous absorption
!                   is represented by a constant value of single scattering
!                   albedo. The default value and the v197 setting are
!                   both false.
!
!      L_anom_abs_v Logical namelist variable. If true, anomalous absorption
!                   is computed as a piecewise continuous function of zonal
!                   mean saturation water vapor mixing ratio at the model's
!                   vertical level closest to the earth's surface. The
!                   default value is false. The analogous namelist variable
!                   in v197, LWALBV is set to TRUE. 
!                   If both L_anom_abs_g and L_anom_abs_v are set to FALSE,
!                   then the single scattering albedo, w0,  assumes its normal
!                   constant value of 0.9942, in tg's version. w0 is a
!                   variable in Steve Klein's version.
!                   L_anom_abs_v takes precedence over L_anom_abs_g, if
!                   the former is TRUE and the latter is FALSE.
!
!     qmix_kx    Zonal mean saturation water vapor mixing ratio at the
!                   model's vertical level closest to the earth's surface.
!                   It is more convenient to pass the zonal mean, in case
!                   single column tests are performed. Alternatively,
!                   the zonal mean could be computed within this subroutine. 
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!       ------------
!INPUT/OUTPUT:
!       ------------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emissivity

!            NOTE:  In tg's version, LWP and IWP are effective cloud
!                   water paths. They could be computed either in this
!                   subroutine or in subroutine cloud_water_path.

!      LWP          cloud liquid water path (kg of condensate per square meter)
!      IWP          cloud ice path (kg of condensate per square meter)
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq      optical depth            at each band for cloud liquid
!       tau_ice      optical depth            at each band for cloud ice
!       w0_liq       single scattering albedo at each band for cloud liquid
!       w0_ice       single scattering albedo at each band for cloud ice
!       g_liq        asymmetry parameter      at each band for cloud liquid
!       g_ice        asymmetry parameter      at each band for cloud ice
!
!
!                   In Tony Gordon's v197 version only.
!
!       CWP           total cloud water path, i.e., cloud liquid plus
!                        cloud ice (kg of condensate per square meter)
!
!                    Parameters (defined above)
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!       qsat_min      Minimum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!       qsat_trans    Transition value of zonal mean saturation mixing ratio
!                        between two branches of a piecewise continuous
!                        function entering into the variable anomalous
!                        absorption single scattering albedo calculations.
!       qsat_max      Maximum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!
!       w0_norm_uv    Normal, constant single scattering albedo in the uv-vis
!                        wavelength band of the radiation spectrum.
!       w0_norm_nir   Normal, constant single scattering albedo in the nir
!                        wavelength band of the radiation spectrum.
!       w0_anom1_nir  Asymptotic minimum value of single scattering albedo
!                        for variable anomalous absorption; also the
!                        low constant value, if L_anom_abs_g is set to TRUE.
!       w0_anom2_nir  Asymptotic maximum value of single scattering albedo
!                        for variable anomalous absorption,usually occurring
!                        at high latitudes.
! 
!       g_norm        Normal, constant asymmetry parameter used in Tony Gordon's
!                        v197 scheme. It is independent of wavelength.                                
!

!      MAX_BAND     maximum number of wave number bands for tau,etc.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!  in tony gordon's v197 scheme only:
!
real,     intent (in),     dimension(:,:)      :: qmix_kx
real,     intent (in),     dimension(:,:,:)    :: tempcld


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real,     intent (in   ),dimension(:,:,:,:)   :: tau
real,     intent (out),dimension(:,:,:,:)   :: w0,gg
real,     intent (out),dimension(:,:,:), optional     :: em_lw

real,     intent (out)   ,dimension(:,:,:)  :: LWP,IWP


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  Internal variables
!  ------------------

real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: tau_liq, tau_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: w0_liq, w0_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: g_liq, g_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: wgt_ice, wgt_liq
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: cwp
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: w0_anom_work
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: qmix_kx_work
integer                                                :: k, kmax
integer                                                :: n, max_band

!  Declare tau_chk to compare with tau, if code needs to be debugged.
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: tau_chk

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!
! Code
! ----

!  define maximum number of wave number bands for tau
      max_band = size(tau,4)

!  For tg v197, tau is previously computed in subroutine cloud_optical_depths,
!  while w0 and gg will be reinitialized.
!  Also, internal variables w0_liq, w0_ice, g_liq, and g_ice
!  will be re-initialized, while tau_liq and tau_ice will not be.

!  Therefore, the only output variable to be initialized is em_lw.
!  Comment out the other reinitialization commands.


!  These are Tony Gordon's reinitialized values.
        gg(:,:,:,:)    = 0.85
        w0(:,:,:,1)    = 0.99999
        w0(:,:,:,2:4)  = 0.9942
if (present (em_lw)) then
        em_lw(:,:,:)   = 0.
endif

        w0_liq(:,:,:,1)   = 0.99999
        w0_liq(:,:,:,2:4) = 0.9942
        w0_ice(:,:,:,1)   = 0.99999
        w0_ice(:,:,:,2:4) = 0.9942
        g_liq(:,:,:,:)    = 0.85
        g_ice(:,:,:,:)    = 0.85

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!       Comment out Steve Klein's reinitialized values.
!       tau(:,:,:,:) = 0.
!       gg(:,:,:,:)  = 0.85
!       w0(:,:,:,:)  = 0.95
!       em_lw(:,:,:) = 0.

!       w0_liq(:,:,:,:) = 0.95
!       w0_ice(:,:,:,:) = 0.95
!       tau_liq(:,:,:,:)= 0.
!       tau_ice(:,:,:,:)= 0.



    !---------------   COMPUTE OPTICAL DEPTH ---------------------------!



        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately


!       by Tony Gordon's v197 scheme

        WHERE (tempcld(:,:,:) .le. tk_all_ice)
               wgt_liq(:,:,:) = 0.
               wgt_ice(:,:,:) = 1.
        END WHERE

        WHERE (tempcld(:,:,:) .ge. tk_all_liq)
               wgt_liq(:,:,:) = 1.
               wgt_ice(:,:,:) = 0.
        END WHERE

        WHERE (tempcld(:,:,:) .gt. tk_all_ice .and. tempcld(:,:,:) &
               .lt. tk_all_liq)
               wgt_liq(:,:,:) = (tempcld(:,:,:) - tk_all_ice) / &
                                (tk_all_liq - tk_all_ice)
               wgt_ice(:,:,:) = 1. - wgt_liq(:,:,:)
        END WHERE
                     
        CWP(:,:,:) = tau(:,:,:,1) / &
                     (k_sw_liq * wgt_liq(:,:,:) + &
                      k_sw_ice * wgt_ice(:,:,:) )

        LWP(:,:,:) = wgt_liq(:,:,:) * CWP(:,:,:)
        IWP(:,:,:) = wgt_ice(:,:,:) * CWP(:,:,:)

        tau_liq(:,:,:,1)   = k_sw_liq * LWP(:,:,:)
        tau_ice(:,:,:,1)   = k_sw_ice * IWP(:,:,:)

!  tau_liq and tau_ice are purely diagnostic, since tau is already known.
!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes.

             if (max_band .ge. 2) then
        do n=2,max_band
        tau_liq(:,:,:,n) = tau_liq(:,:,:,1)
        tau_ice(:,:,:,n) = tau_ice(:,:,:,1)
        end do
             endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Compute total cloud optical depth, using same formula as Steve Klein.
!  Note:  Comment out the following command in Tony Gordon's v197 scheme.
!         tau should have the same as the input, except for roundoff error.

!         tau(:,:,:,:)     = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

!  Define tau_chk to compare with tau, if code needs to be debugged.
          tau_chk(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



   !---------------   COMPUTE SINGLE SCATTERING ALBEDO ----------------!


!
        w0_liq(:,:,:,1) =  w0_norm_uv
        w0_ice(:,:,:,1) =  w0_norm_uv

     IF (.not. L_anom_abs_g .and. .not. L_anom_abs_v) THEN

!       Specify tg's normal single scattering albedos for NIR bands.

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_norm_nir
        w0_ice(:,:,:,n) =  w0_norm_nir
        end do
             endif

     ENDIF

     IF (L_anom_abs_g) THEN

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_anom1_nir
        w0_ice(:,:,:,n) =  w0_anom1_nir
        end do
             endif

     ENDIF

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!    Broadcast qmix_kx to generate WHERE statements with correct syntax.

     kmax = SIZE(LWP,3)

     DO k = 1,kmax
        qmix_kx_work(:,:,k) = qmix_kx(:,:)
     END DO
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! 
     IF (L_anom_abs_v) THEN

       WHERE (qmix_kx_work(:,:,:) .le. qsat_min)

!      Apply lower asymptotic limit to anomalously weak cloud absorption,
!      i.e., upper asymptotic limit of w0.

!      This situation should not occur for saturation mixing ratios. However,
!      negative mixing ratios are possible in the spectral AGCM,
!      especially in cold temperature, e.g., high latitude regions.
!      Retain this WHERE loop, in case parameterization is ever changed 
!      from saturation mixing ratio to mixing ratio.
      
              w0_anom_work(:,:,:) = w0_anom2_nir

       END WHERE 

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_min .and.                  &
              qmix_kx_work(:,:,:) .lt. qsat_trans)

!      Anomalously weak cloud absorption relative to the reference value 
!      w0_norm_nir will tend to occur at higher latitudes for these values
!      of w0.

            w0_anom_work(:,:,:) = w0_anom2_nir -                    &
               (w0_anom2_nir - w0_norm_nir)    *                    &
               ( (qmix_kx_work(:,:,:) - qsat_min) / (qsat_trans - qsat_min))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .eq. qsat_trans)

!      The reference value of nir single scattering albedo will be used.

           w0_anom_work(:,:,:) = w0_norm_nir

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_trans .and.                &
              qmix_kx_work(:,:,:) .lt. qsat_max)

!      Anomalously high absorption relative to the reference value w0_norm_nir
!      will tend to occur at tropical and subtropical latitudes for 
!      these values of w0.

           w0_anom_work(:,:,:) = w0_norm_nir  -                    &
              (w0_norm_nir - w0_anom1_nir)    *                    &
              ( (qmix_kx_work(:,:,:) - qsat_trans) / (qsat_max - qsat_trans))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .ge. qsat_max)

!      Apply upper asymptotic limit to the anomalous absorption, i.e.,
!      lower asymptotic limit of w0.

              w0_anom_work = w0_anom1_nir

       END WHERE

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
       w0_liq(:,:,:,n) = w0_anom_work(:,:,:)
       w0_ice(:,:,:,n) = w0_anom_work(:,:,:)
        end do
             endif

     ENDIF



! compute total single scattering albedo
        WHERE (tau(:,:,:,:) .gt. 0.)
               w0(:,:,:,:) = ( w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                               w0_ice(:,:,:,:) * tau_ice(:,:,:,:) ) / &
                             tau(:,:,:,:)
        END WHERE

   !---------------   COMPUTE ASYMMETRY PARAMETER --------------------!

        WHERE (tau(:,:,:,:) .gt. 0. )
              gg(:,:,:,:) = ( &
                 w0_liq(:,:,:,:) * g_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                 w0_ice(:,:,:,:) * g_ice(:,:,:,:) * tau_ice(:,:,:,:) ) &
                       /          (w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                                   w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )
        END WHERE


   !---------------   COMPUTE LONGWAVE EMISSIVITY --------------------!


       if (present(em_lw)) then
!
!  In Tony Gordon's v197 scheme, k_lw-liq and k_lw_ice are parameters.
!        k_lw_liq(:,:,:) = 140.
!        k_lw_ice(:,:,:) = 100.

! compute combined emmisivity
        em_lw(:,:,:) =  1. - exp( -1. * ( k_lw_liq * LWP(:,:,:) + &
                                          k_lw_ice * IWP(:,:,:) ) )

     endif

   !--------------    RANGE LIMIT QUANTITIES --------------------------!

!       WHERE (tau(:,:,:,:) .lt. taumin)
!              tau(:,:,:,:) = taumin
!       END WHERE


end subroutine CLOUD_OPT_PROP_tg


subroutine CLOUD_OPT_PROP_tg3(LWP,IWP,       &
                        tau,w0,gg,          &
                        tempcld,qmix_kx, em_lw )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau    :optical depth in each band
!               2. w0     :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw  :longwave cloud emissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
! *************************    WARNING    *****************************
!
!   The above bands are used by Steve Klein.
!   We retain the scheme from the v197 frozen model,instead.
!   Nominally, our band 2 is expanded into bands 2 + 3 + 4 of Slingo.
!   The same cloud optical depth is specified in all 4 bands.
!   The same asymmetry parameter is specified in all 4 bands.
!   For single scattering albedo w0, the uv value is specified in band 1,
!   while the nir value is specified in bands 2, 3, and 4.

!   ****  WARNING    ****  The code is intended to be applied to 2 to 4 bands.
!   An error check to check that this condition is satisfied is advised.  

! *********************************************************************
!            BAND               v197     
!             1               0.25-0.70
!             2               0.70-4.00
! *********************************************************************
!
!

!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = 0.9942  (in v197 - standard)
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = F(Z.M. mixing ratio at lowest model level)
!                                    (in v197 - anomalous absorption)
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!           g(:,:,1:2) = 0.85 in v197
!   
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation could be rewritten, after algebraic manipulation, as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    However, the other form of the equation, i.e., 
!    1 - exp(tau_liq + tau_ice) will actually be solved.

! *******************************************************************
!
!
!   (6)  v197 only: Must first solve for LWP and IWP knowing
!                   tau, k_sw_liq, k_sw_ice, wgt_liq and wgt_ice.
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!


!      L_anom_abs-g Logical namelist variable. If true, anomalous absorption
!                   is represented by a constant value of single scattering
!                   albedo. The default value and the v197 setting are
!                   both false.
!
!      L_anom_abs_v Logical namelist variable. If true, anomalous absorption
!                   is computed as a piecewise continuous function of zonal
!                   mean saturation water vapor mixing ratio at the model's
!                   vertical level closest to the earth's surface. The
!                   default value is false. The analogous namelist variable
!                   in v197, LWALBV is set to TRUE. 
!                   If both L_anom_abs_g and L_anom_abs_v are set to FALSE,
!                   then the single scattering albedo, w0,  assumes its normal
!                   constant value of 0.9942, in tg's version. w0 is a
!                   variable in Steve Klein's version.
!                   L_anom_abs_v takes precedence over L_anom_abs_g, if
!                   the former is TRUE and the latter is FALSE.
!
!     qmix_kx    Zonal mean saturation water vapor mixing ratio at the
!                   model's vertical level closest to the earth's surface.
!                   It is more convenient to pass the zonal mean, in case
!                   single column tests are performed. Alternatively,
!                   the zonal mean could be computed within this subroutine. 
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!       ------------
!INPUT/OUTPUT:
!       ------------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emissivity

!            NOTE:  In tg's version, LWP and IWP are effective cloud
!                   water paths. They could be computed either in this
!                   subroutine or in subroutine cloud_water_path.

!      LWP          cloud liquid water path (kg of condensate per square meter)
!      IWP          cloud ice path (kg of condensate per square meter)
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq      optical depth            at each band for cloud liquid
!       tau_ice      optical depth            at each band for cloud ice
!       w0_liq       single scattering albedo at each band for cloud liquid
!       w0_ice       single scattering albedo at each band for cloud ice
!       g_liq        asymmetry parameter      at each band for cloud liquid
!       g_ice        asymmetry parameter      at each band for cloud ice
!
!
!                   In Tony Gordon's v197 version only.
!
!       CWP           total cloud water path, i.e., cloud liquid plus
!                        cloud ice (kg of condensate per square meter)
!
!                    Parameters (defined above)
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!       qsat_min      Minimum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!       qsat_trans    Transition value of zonal mean saturation mixing ratio
!                        between two branches of a piecewise continuous
!                        function entering into the variable anomalous
!                        absorption single scattering albedo calculations.
!       qsat_max      Maximum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!
!       w0_norm_uv    Normal, constant single scattering albedo in the uv-vis
!                        wavelength band of the radiation spectrum.
!       w0_norm_nir   Normal, constant single scattering albedo in the nir
!                        wavelength band of the radiation spectrum.
!       w0_anom1_nir  Asymptotic minimum value of single scattering albedo
!                        for variable anomalous absorption; also the
!                        low constant value, if L_anom_abs_g is set to TRUE.
!       w0_anom2_nir  Asymptotic maximum value of single scattering albedo
!                        for variable anomalous absorption,usually occurring
!                        at high latitudes.
! 
!       g_norm        Normal, constant asymmetry parameter used in Tony Gordon's
!                        v197 scheme. It is independent of wavelength.                                
!

!      MAX_BAND     maximum number of wave number bands for tau,etc.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!  in tony gordon's v197 scheme only:
!
real,     intent (in),     dimension(:,:)      :: qmix_kx
real,     intent (in),     dimension(:,:,:)    :: tempcld


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real,     intent (in   ),dimension(:,:,:,:)   :: tau
real,     intent (out),dimension(:,:,:,:)   :: w0,gg
real,     intent (out),dimension(:,:,:), optional     :: em_lw

real,     intent ( in)   ,dimension(:,:,:)  :: LWP,IWP


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  Internal variables
!  ------------------

real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: tau_liq, tau_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: w0_liq, w0_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: g_liq, g_ice
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: w0_anom_work
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3))   :: qmix_kx_work
integer                                                :: k, kmax
integer                                                :: n, max_band

!  Declare tau_chk to compare with tau, if code needs to be debugged.
real, dimension(size(lwp,1),size(lwp,2),size(lwp,3),4) :: tau_chk

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!
! Code
! ----

!  define maximum number of wave number bands for tau
      max_band = size(tau,4)

!  For tg v197, tau is previously computed in subroutine cloud_optical_depths,
!  while w0 and gg will be reinitialized.
!  Also, internal variables w0_liq, w0_ice, g_liq, and g_ice
!  will be re-initialized, while tau_liq and tau_ice will not be.

!  Therefore, the only output variable to be initialized is em_lw.
!  Comment out the other reinitialization commands.


!  These are Tony Gordon's reinitialized values.
        gg(:,:,:,:)    = 0.85
        w0(:,:,:,1)    = 0.99999
        w0(:,:,:,2:4)  = 0.9942
if (present (em_lw)) then
        em_lw(:,:,:)   = 0.
endif

        w0_liq(:,:,:,1)   = 0.99999
        w0_liq(:,:,:,2:4) = 0.9942
        w0_ice(:,:,:,1)   = 0.99999
        w0_ice(:,:,:,2:4) = 0.9942
        g_liq(:,:,:,:)    = 0.85
        g_ice(:,:,:,:)    = 0.85

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!       Comment out Steve Klein's reinitialized values.
!       tau(:,:,:,:) = 0.
!       gg(:,:,:,:)  = 0.85
!       w0(:,:,:,:)  = 0.95
!       em_lw(:,:,:) = 0.

!       w0_liq(:,:,:,:) = 0.95
!       w0_ice(:,:,:,:) = 0.95
!       tau_liq(:,:,:,:)= 0.
!       tau_ice(:,:,:,:)= 0.



    !---------------   COMPUTE OPTICAL DEPTH ---------------------------!



        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately


!       by Tony Gordon's v197 scheme


        tau_liq(:,:,:,1)   = k_sw_liq * LWP(:,:,:)
        tau_ice(:,:,:,1)   = k_sw_ice * IWP(:,:,:)

!  tau_liq and tau_ice are purely diagnostic, since tau is already known.
!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes.

             if (max_band .ge. 2) then
        do n=2,max_band
        tau_liq(:,:,:,n) = tau_liq(:,:,:,1)
        tau_ice(:,:,:,n) = tau_ice(:,:,:,1)
        end do
             endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Compute total cloud optical depth, using same formula as Steve Klein.
!  Note:  Comment out the following command in Tony Gordon's v197 scheme.
!         tau should have the same as the input, except for roundoff error.

!         tau(:,:,:,:)     = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

!  Define tau_chk to compare with tau, if code needs to be debugged.
          tau_chk(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



   !---------------   COMPUTE SINGLE SCATTERING ALBEDO ----------------!


!
        w0_liq(:,:,:,1) =  w0_norm_uv
        w0_ice(:,:,:,1) =  w0_norm_uv

     IF (.not. L_anom_abs_g .and. .not. L_anom_abs_v) THEN

!       Specify tg's normal single scattering albedos for NIR bands.

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_norm_nir
        w0_ice(:,:,:,n) =  w0_norm_nir
        end do
             endif

     ENDIF

     IF (L_anom_abs_g) THEN

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_anom1_nir
        w0_ice(:,:,:,n) =  w0_anom1_nir
        end do
             endif

     ENDIF

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!    Broadcast qmix_kx to generate WHERE statements with correct syntax.

     kmax = SIZE(LWP,3)

     DO k = 1,kmax
        qmix_kx_work(:,:,k) = qmix_kx(:,:)
     END DO
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! 
     IF (L_anom_abs_v) THEN

       WHERE (qmix_kx_work(:,:,:) .le. qsat_min)

!      Apply lower asymptotic limit to anomalously weak cloud absorption,
!      i.e., upper asymptotic limit of w0.

!      This situation should not occur for saturation mixing ratios. However,
!      negative mixing ratios are possible in the spectral AGCM,
!      especially in cold temperature, e.g., high latitude regions.
!      Retain this WHERE loop, in case parameterization is ever changed 
!      from saturation mixing ratio to mixing ratio.
      
              w0_anom_work(:,:,:) = w0_anom2_nir

       END WHERE 

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_min .and.                  &
              qmix_kx_work(:,:,:) .lt. qsat_trans)

!      Anomalously weak cloud absorption relative to the reference value 
!      w0_norm_nir will tend to occur at higher latitudes for these values
!      of w0.

            w0_anom_work(:,:,:) = w0_anom2_nir -                    &
               (w0_anom2_nir - w0_norm_nir)    *                    &
               ( (qmix_kx_work(:,:,:) - qsat_min) / (qsat_trans - qsat_min))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .eq. qsat_trans)

!      The reference value of nir single scattering albedo will be used.

           w0_anom_work(:,:,:) = w0_norm_nir

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_trans .and.                &
              qmix_kx_work(:,:,:) .lt. qsat_max)

!      Anomalously high absorption relative to the reference value w0_norm_nir
!      will tend to occur at tropical and subtropical latitudes for 
!      these values of w0.

           w0_anom_work(:,:,:) = w0_norm_nir  -                    &
              (w0_norm_nir - w0_anom1_nir)    *                    &
              ( (qmix_kx_work(:,:,:) - qsat_trans) / (qsat_max - qsat_trans))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .ge. qsat_max)

!      Apply upper asymptotic limit to the anomalous absorption, i.e.,
!      lower asymptotic limit of w0.

              w0_anom_work = w0_anom1_nir

       END WHERE

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
       w0_liq(:,:,:,n) = w0_anom_work(:,:,:)
       w0_ice(:,:,:,n) = w0_anom_work(:,:,:)
        end do
             endif

     ENDIF



! compute total single scattering albedo
        WHERE (tau(:,:,:,:) .gt. 0.)
               w0(:,:,:,:) = ( w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                               w0_ice(:,:,:,:) * tau_ice(:,:,:,:) ) / &
                             tau(:,:,:,:)
        END WHERE

   !---------------   COMPUTE ASYMMETRY PARAMETER --------------------!

        WHERE (tau(:,:,:,:) .gt. 0. )
              gg(:,:,:,:) = ( &
                 w0_liq(:,:,:,:) * g_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                 w0_ice(:,:,:,:) * g_ice(:,:,:,:) * tau_ice(:,:,:,:) ) &
                       /          (w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                                   w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )
        END WHERE


   !---------------   COMPUTE LONGWAVE EMISSIVITY --------------------!


       if (present(em_lw)) then
!
!  In Tony Gordon's v197 scheme, k_lw-liq and k_lw_ice are parameters.
!        k_lw_liq(:,:,:) = 140.
!        k_lw_ice(:,:,:) = 100.

! compute combined emmisivity
        em_lw(:,:,:) =  1. - exp( -1. * ( k_lw_liq * LWP(:,:,:) + &
                                          k_lw_ice * IWP(:,:,:) ) )

     endif

   !--------------    RANGE LIMIT QUANTITIES --------------------------!

!       WHERE (tau(:,:,:,:) .lt. taumin)
!              tau(:,:,:,:) = taumin
!       END WHERE


end subroutine CLOUD_OPT_PROP_tg3

subroutine CLOUD_OPT_PROP_tg_lw(               tau, liq_frac, em_lw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau    :optical depth in each band
!               2. w0     :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw  :longwave cloud emissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
! *************************    WARNING    *****************************
!
!   The above bands are used by Steve Klein.
!   We retain the scheme from the v197 frozen model,instead.
!   Nominally, our band 2 is expanded into bands 2 + 3 + 4 of Slingo.
!   The same cloud optical depth is specified in all 4 bands.
!   The same asymmetry parameter is specified in all 4 bands.
!   For single scattering albedo w0, the uv value is specified in band 1,
!   while the nir value is specified in bands 2, 3, and 4.

!   ****  WARNING    ****  The code is intended to be applied to 2 to 4 bands.
!   An error check to check that this condition is satisfied is advised.  

! *********************************************************************
!            BAND               v197     
!             1               0.25-0.70
!             2               0.70-4.00
! *********************************************************************
!
!

!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = 0.9942  (in v197 - standard)
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = F(Z.M. mixing ratio at lowest model level)
!                                    (in v197 - anomalous absorption)
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!           g(:,:,1:2) = 0.85 in v197
!   
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation could be rewritten, after algebraic manipulation, as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    However, the other form of the equation, i.e., 
!    1 - exp(tau_liq + tau_ice) will actually be solved.

! *******************************************************************
!
!
!   (6)  v197 only: Must first solve for LWP and IWP knowing
!                   tau, k_sw_liq, k_sw_ice, wgt_liq and wgt_ice.
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!


!      L_anom_abs-g Logical namelist variable. If true, anomalous absorption
!                   is represented by a constant value of single scattering
!                   albedo. The default value and the v197 setting are
!                   both false.
!
!      L_anom_abs_v Logical namelist variable. If true, anomalous absorption
!                   is computed as a piecewise continuous function of zonal
!                   mean saturation water vapor mixing ratio at the model's
!                   vertical level closest to the earth's surface. The
!                   default value is false. The analogous namelist variable
!                   in v197, LWALBV is set to TRUE. 
!                   If both L_anom_abs_g and L_anom_abs_v are set to FALSE,
!                   then the single scattering albedo, w0,  assumes its normal
!                   constant value of 0.9942, in tg's version. w0 is a
!                   variable in Steve Klein's version.
!                   L_anom_abs_v takes precedence over L_anom_abs_g, if
!                   the former is TRUE and the latter is FALSE.
!
!     qmix_kx    Zonal mean saturation water vapor mixing ratio at the
!                   model's vertical level closest to the earth's surface.
!                   It is more convenient to pass the zonal mean, in case
!                   single column tests are performed. Alternatively,
!                   the zonal mean could be computed within this subroutine. 
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!       ------------
!INPUT/OUTPUT:
!       ------------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emissivity

!            NOTE:  In tg's version, LWP and IWP are effective cloud
!                   water paths. They could be computed either in this
!                   subroutine or in subroutine cloud_water_path.

!      LWP          cloud liquid water path (kg of condensate per square meter)
!      IWP          cloud ice path (kg of condensate per square meter)
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq      optical depth            at each band for cloud liquid
!       tau_ice      optical depth            at each band for cloud ice
!       w0_liq       single scattering albedo at each band for cloud liquid
!       w0_ice       single scattering albedo at each band for cloud ice
!       g_liq        asymmetry parameter      at each band for cloud liquid
!       g_ice        asymmetry parameter      at each band for cloud ice
!
!
!                   In Tony Gordon's v197 version only.
!
!       CWP           total cloud water path, i.e., cloud liquid plus
!                        cloud ice (kg of condensate per square meter)
!
!                    Parameters (defined above)
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!       qsat_min      Minimum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!       qsat_trans    Transition value of zonal mean saturation mixing ratio
!                        between two branches of a piecewise continuous
!                        function entering into the variable anomalous
!                        absorption single scattering albedo calculations.
!       qsat_max      Maximum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!
!       w0_norm_uv    Normal, constant single scattering albedo in the uv-vis
!                        wavelength band of the radiation spectrum.
!       w0_norm_nir   Normal, constant single scattering albedo in the nir
!                        wavelength band of the radiation spectrum.
!       w0_anom1_nir  Asymptotic minimum value of single scattering albedo
!                        for variable anomalous absorption; also the
!                        low constant value, if L_anom_abs_g is set to TRUE.
!       w0_anom2_nir  Asymptotic maximum value of single scattering albedo
!                        for variable anomalous absorption,usually occurring
!                        at high latitudes.
! 
!       g_norm        Normal, constant asymmetry parameter used in Tony Gordon's
!                        v197 scheme. It is independent of wavelength.                                
!

!      MAX_BAND     maximum number of wave number bands for tau,etc.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!  in tony gordon's v197 scheme only:
!
!real,     intent (in),     dimension(:,:)      :: qmix_kx
!real,     intent (in),     dimension(:,:,:)    :: tempcld


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real,     intent (out),dimension(:,:,:)     :: em_lw
real,     intent ( in)   ,dimension(:,:,:)  ::                liq_frac
real,     intent ( in)   ,dimension(:,:,:,:)  :: tau                 


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  Internal variables
!  ------------------

real, dimension(size(tau   ,1),size(tau   ,2),size(tau   ,3))   :: cwp, lwp, iwp

!  Declare tau_chk to compare with tau, if code needs to be debugged.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!
! Code
! ----

!  define maximum number of wave number bands for tau
!      max_band = size(tau,4)

!  For tg v197, tau is previously computed in subroutine cloud_optical_depths,
!  while w0 and gg will be reinitialized.
!  Also, internal variables w0_liq, w0_ice, g_liq, and g_ice
!  will be re-initialized, while tau_liq and tau_ice will not be.

!  Therefore, the only output variable to be initialized is em_lw.
!  Comment out the other reinitialization commands.


        em_lw(:,:,:)   = 0.


! compute combined emmisivity
        CWP(:,:,:) = tau(:,:,:,1) / &
                     (k_sw_liq * liq_frac(:,:,:) + &
                      k_sw_ice * (1.0-liq_frac(:,:,:)) )

        LWP(:,:,:) = liq_frac(:,:,:) * CWP(:,:,:)
        IWP(:,:,:) = (1.0-liq_frac(:,:,:)) * CWP(:,:,:)
        em_lw(:,:,:) =  1. - exp( -1. * ( k_lw_liq * LWP(:,:,:) + &
                                          k_lw_ice * IWP(:,:,:) ) )


   !--------------    RANGE LIMIT QUANTITIES --------------------------!



end subroutine CLOUD_OPT_PROP_tg_lw

subroutine CLOUD_OPT_PROP_tg_sw(         liq_frac,       &
                        tau, direct,              &
                                qmix_kx, cosz, cuvrf, cirrf, &
                        cuvab, cirab)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau    :optical depth in each band
!               2. w0     :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw  :longwave cloud emissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
! *************************    WARNING    *****************************
!
!   The above bands are used by Steve Klein.
!   We retain the scheme from the v197 frozen model,instead.
!   Nominally, our band 2 is expanded into bands 2 + 3 + 4 of Slingo.
!   The same cloud optical depth is specified in all 4 bands.
!   The same asymmetry parameter is specified in all 4 bands.
!   For single scattering albedo w0, the uv value is specified in band 1,
!   while the nir value is specified in bands 2, 3, and 4.

!   ****  WARNING    ****  The code is intended to be applied to 2 to 4 bands.
!   An error check to check that this condition is satisfied is advised.  

! *********************************************************************
!            BAND               v197     
!             1               0.25-0.70
!             2               0.70-4.00
! *********************************************************************
!
!

!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = 0.9942  (in v197 - standard)
!   w0(:,:,1) = 0.99999;  w0(:,:,2) = F(Z.M. mixing ratio at lowest model level)
!                                    (in v197 - anomalous absorption)
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!           g(:,:,1:2) = 0.85 in v197
!   
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation could be rewritten, after algebraic manipulation, as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    However, the other form of the equation, i.e., 
!    1 - exp(tau_liq + tau_ice) will actually be solved.

! *******************************************************************
!
!
!   (6)  v197 only: Must first solve for LWP and IWP knowing
!                   tau, k_sw_liq, k_sw_ice, wgt_liq and wgt_ice.
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!


!      L_anom_abs-g Logical namelist variable. If true, anomalous absorption
!                   is represented by a constant value of single scattering
!                   albedo. The default value and the v197 setting are
!                   both false.
!
!      L_anom_abs_v Logical namelist variable. If true, anomalous absorption
!                   is computed as a piecewise continuous function of zonal
!                   mean saturation water vapor mixing ratio at the model's
!                   vertical level closest to the earth's surface. The
!                   default value is false. The analogous namelist variable
!                   in v197, LWALBV is set to TRUE. 
!                   If both L_anom_abs_g and L_anom_abs_v are set to FALSE,
!                   then the single scattering albedo, w0,  assumes its normal
!                   constant value of 0.9942, in tg's version. w0 is a
!                   variable in Steve Klein's version.
!                   L_anom_abs_v takes precedence over L_anom_abs_g, if
!                   the former is TRUE and the latter is FALSE.
!
!     qmix_kx    Zonal mean saturation water vapor mixing ratio at the
!                   model's vertical level closest to the earth's surface.
!                   It is more convenient to pass the zonal mean, in case
!                   single column tests are performed. Alternatively,
!                   the zonal mean could be computed within this subroutine. 
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!       ------------
!INPUT/OUTPUT:
!       ------------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emissivity

!            NOTE:  In tg's version, LWP and IWP are effective cloud
!                   water paths. They could be computed either in this
!                   subroutine or in subroutine cloud_water_path.

!      LWP          cloud liquid water path (kg of condensate per square meter)
!      IWP          cloud ice path (kg of condensate per square meter)
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq      optical depth            at each band for cloud liquid
!       tau_ice      optical depth            at each band for cloud ice
!       w0_liq       single scattering albedo at each band for cloud liquid
!       w0_ice       single scattering albedo at each band for cloud ice
!       g_liq        asymmetry parameter      at each band for cloud liquid
!       g_ice        asymmetry parameter      at each band for cloud ice
!
!
!                   In Tony Gordon's v197 version only.
!
!       CWP           total cloud water path, i.e., cloud liquid plus
!                        cloud ice (kg of condensate per square meter)
!
!                    Parameters (defined above)
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!       qsat_min      Minimum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!       qsat_trans    Transition value of zonal mean saturation mixing ratio
!                        between two branches of a piecewise continuous
!                        function entering into the variable anomalous
!                        absorption single scattering albedo calculations.
!       qsat_max      Maximum value of zonal mean saturation water vapor mixing
!                        ratio entering into the variable anomalous absorption
!                        single scattering albedo calculations.
!
!       w0_norm_uv    Normal, constant single scattering albedo in the uv-vis
!                        wavelength band of the radiation spectrum.
!       w0_norm_nir   Normal, constant single scattering albedo in the nir
!                        wavelength band of the radiation spectrum.
!       w0_anom1_nir  Asymptotic minimum value of single scattering albedo
!                        for variable anomalous absorption; also the
!                        low constant value, if L_anom_abs_g is set to TRUE.
!       w0_anom2_nir  Asymptotic maximum value of single scattering albedo
!                        for variable anomalous absorption,usually occurring
!                        at high latitudes.
! 
!       g_norm        Normal, constant asymmetry parameter used in Tony Gordon's
!                        v197 scheme. It is independent of wavelength.                                
!

!      MAX_BAND     maximum number of wave number bands for tau,etc.

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!  in tony gordon's v197 scheme only:
!
real,     intent (in),     dimension(:,:)      :: qmix_kx, cosz
logical,  intent (in),     dimension(:,:,:)    :: direct           


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real,     intent (in),dimension(:,:,:,:)   :: tau

real,     intent ( in)   ,dimension(:,:,:)  ::          liq_frac
real,     intent (out),     dimension(:,:,:)    :: cuvrf, cirrf,   &
                                                   cuvab, cirab


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  Internal variables
!  ------------------

real, dimension(size(tau,1),size(tau,2),size(tau,3),4) :: tau_liq, tau_ice
real, dimension(size(tau,1),size(tau,2),size(tau,3),4) :: w0_liq, w0_ice
real, dimension(size(tau,1),size(tau,2),size(tau,3),4) :: g_liq, g_ice
real, dimension(size(tau,1),size(tau,2),size(tau,3))   :: w0_anom_work
real, dimension(size(tau,1),size(tau,2),size(tau,3))   :: qmix_kx_work
real, dimension(size(tau,1),size(tau,2),size(tau,3))   :: lwp_new, iwp_new, cwp
integer                                                :: k, kmax
integer                                                :: n, max_band

!  Declare tau_chk to compare with tau, if code needs to be debugged.
real, dimension(size(tau,1),size(tau,2),size(tau,3),4) :: tau_chk
real, dimension(size(tau,1),size(tau,2),size(tau,3),   &
                size(tau,4)) :: w0, gg 

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!
! Code
! ----

!  define maximum number of wave number bands for tau
      max_band = size(tau,4)

!  For tg v197, tau is previously computed in subroutine cloud_optical_depths,
!  while w0 and gg will be reinitialized.
!  Also, internal variables w0_liq, w0_ice, g_liq, and g_ice
!  will be re-initialized, while tau_liq and tau_ice will not be.

!  Therefore, the only output variable to be initialized is em_lw.
!  Comment out the other reinitialization commands.


!  These are Tony Gordon's reinitialized values.
        gg(:,:,:,:)    = 0.85
        w0(:,:,:,1)    = 0.99999
        w0(:,:,:,2:4)  = 0.9942

        w0_liq(:,:,:,1)   = 0.99999
        w0_liq(:,:,:,2:4) = 0.9942
        w0_ice(:,:,:,1)   = 0.99999
        w0_ice(:,:,:,2:4) = 0.9942
        g_liq(:,:,:,:)    = 0.85
        g_ice(:,:,:,:)    = 0.85

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!       Comment out Steve Klein's reinitialized values.
!       tau(:,:,:,:) = 0.
!       gg(:,:,:,:)  = 0.85
!       w0(:,:,:,:)  = 0.95
!       em_lw(:,:,:) = 0.

!       w0_liq(:,:,:,:) = 0.95
!       w0_ice(:,:,:,:) = 0.95
!       tau_liq(:,:,:,:)= 0.
!       tau_ice(:,:,:,:)= 0.



    !---------------   COMPUTE OPTICAL DEPTH ---------------------------!



        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately


                     
        CWP(:,:,:) = tau(:,:,:,1) / &
                     (k_sw_liq * liq_frac(:,:,:) + &
                      k_sw_ice * (1.-liq_frac(:,:,:)) )

        LWP_new(:,:,:) = liq_frac(:,:,:) * CWP(:,:,:)
        IWP_new(:,:,:) = (1.0-liq_frac(:,:,:)) * CWP(:,:,:)

        tau_liq(:,:,:,1)   = k_sw_liq * LWP_new(:,:,:)
        tau_ice(:,:,:,1)   = k_sw_ice * IWP_new(:,:,:)

!  tau_liq and tau_ice are purely diagnostic, since tau is already known.
!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes.

             if (max_band .ge. 2) then
        do n=2,max_band
        tau_liq(:,:,:,n) = tau_liq(:,:,:,1)
        tau_ice(:,:,:,n) = tau_ice(:,:,:,1)
        end do
             endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Compute total cloud optical depth, using same formula as Steve Klein.
!  Note:  Comment out the following command in Tony Gordon's v197 scheme.
!         tau should have the same as the input, except for roundoff error.

!         tau(:,:,:,:)     = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

!  Define tau_chk to compare with tau, if code needs to be debugged.
          tau_chk(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



   !---------------   COMPUTE SINGLE SCATTERING ALBEDO ----------------!


!
        w0_liq(:,:,:,1) =  w0_norm_uv
        w0_ice(:,:,:,1) =  w0_norm_uv

     IF (.not. L_anom_abs_g .and. .not. L_anom_abs_v) THEN

!       Specify tg's normal single scattering albedos for NIR bands.

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_norm_nir
        w0_ice(:,:,:,n) =  w0_norm_nir
        end do
             endif

     ENDIF

     IF (L_anom_abs_g) THEN

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
        w0_liq(:,:,:,n) =  w0_anom1_nir
        w0_ice(:,:,:,n) =  w0_anom1_nir
        end do
             endif

     ENDIF

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!    Broadcast qmix_kx to generate WHERE statements with correct syntax.

     kmax = SIZE(tau,3)

     DO k = 1,kmax
        qmix_kx_work(:,:,k) = qmix_kx(:,:)
     END DO
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! 
     IF (L_anom_abs_v) THEN

       WHERE (qmix_kx_work(:,:,:) .le. qsat_min)

!      Apply lower asymptotic limit to anomalously weak cloud absorption,
!      i.e., upper asymptotic limit of w0.

!      This situation should not occur for saturation mixing ratios. However,
!      negative mixing ratios are possible in the spectral AGCM,
!      especially in cold temperature, e.g., high latitude regions.
!      Retain this WHERE loop, in case parameterization is ever changed 
!      from saturation mixing ratio to mixing ratio.
      
              w0_anom_work(:,:,:) = w0_anom2_nir

       END WHERE 

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_min .and.                  &
              qmix_kx_work(:,:,:) .lt. qsat_trans)

!      Anomalously weak cloud absorption relative to the reference value 
!      w0_norm_nir will tend to occur at higher latitudes for these values
!      of w0.

            w0_anom_work(:,:,:) = w0_anom2_nir -                    &
               (w0_anom2_nir - w0_norm_nir)    *                    &
               ( (qmix_kx_work(:,:,:) - qsat_min) / (qsat_trans - qsat_min))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .eq. qsat_trans)

!      The reference value of nir single scattering albedo will be used.

           w0_anom_work(:,:,:) = w0_norm_nir

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .gt. qsat_trans .and.                &
              qmix_kx_work(:,:,:) .lt. qsat_max)

!      Anomalously high absorption relative to the reference value w0_norm_nir
!      will tend to occur at tropical and subtropical latitudes for 
!      these values of w0.

           w0_anom_work(:,:,:) = w0_norm_nir  -                    &
              (w0_norm_nir - w0_anom1_nir)    *                    &
              ( (qmix_kx_work(:,:,:) - qsat_trans) / (qsat_max - qsat_trans))

       END WHERE

       WHERE (qmix_kx_work(:,:,:) .ge. qsat_max)

!      Apply upper asymptotic limit to the anomalous absorption, i.e.,
!      lower asymptotic limit of w0.

              w0_anom_work = w0_anom1_nir

       END WHERE

!  Generalize code to n bands, though it may need to be revised,
!  if max_band changes

             if (max_band .ge. 2) then
        do n=2,max_band
       w0_liq(:,:,:,n) = w0_anom_work(:,:,:)
       w0_ice(:,:,:,n) = w0_anom_work(:,:,:)
        end do
             endif

     ENDIF



! compute total single scattering albedo
        WHERE (tau(:,:,:,:) .gt. 0.)
               w0(:,:,:,:) = ( w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                               w0_ice(:,:,:,:) * tau_ice(:,:,:,:) ) / &
                             tau(:,:,:,:)
        END WHERE

   !---------------   COMPUTE ASYMMETRY PARAMETER --------------------!

        WHERE (tau(:,:,:,:) .gt. 0. )
              gg(:,:,:,:) = ( &
                 w0_liq(:,:,:,:) * g_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                 w0_ice(:,:,:,:) * g_ice(:,:,:,:) * tau_ice(:,:,:,:) ) &
                       /          (w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                                   w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )
        END WHERE



       call cloud_rad_k_diag(tau, direct, w0,gg,cosz,cuvrf,cirrf,cuvab,cirab )

end subroutine CLOUD_OPT_PROP_tg_sw




subroutine CLOUD_OPT_PROP_tg2 (tau, tempcld,           liq_frac)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!       ------------
!INPUT/OUTPUT:
!       ------------
!
!      tau          optical depth in each band

!            NOTE:  In tg's version, LWP and IWP are effective cloud
!                   water paths. They could be computed either in this
!                   subroutine or in subroutine cloud_water_path.

!      LWP          cloud liquid water path (kg of condensate per square meter)
!      IWP          cloud ice path (kg of condensate per square meter)
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq      optical depth            at each band for cloud liquid
!       tau_ice      optical depth            at each band for cloud ice
!
!
!                   In Tony Gordon's v197 version only.
!
!       CWP           total cloud water path, i.e., cloud liquid plus
!                        cloud ice (kg of condensate per square meter)
!
!                    Parameters (defined above)
!
!       k_lw_liq     liquid cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_lw_ice     ice cloud mass absorption coefficient for longwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_liq     liquid cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!       k_sw_ice     ice cloud mass absorption coefficient for shortwave
!                       portion of the spectrum (meters**2./kg of condensate)
!
!       tk_all_ice    minimum temperature at which cloud liquid water phase
!                       can exist (degrees Kelvin)
!       tk_all_liq    maximum temperature at which cloud ice can exist
!                       (degrees Kelvin)
!       wgt_liq       The ratio of liquid water path to total cloud water path,
!                        i.e., LWP / CWP
!       wgt_ice       The ratio of ice water path to total cloud water path,
!                        i.e., IWP / CWP
!
!

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real,     intent (in   ),dimension(:,:,:,:)   :: tau
real,     intent (in),     dimension(:,:,:)    :: tempcld


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


real,     intent (out)   ,dimension(:,:,:)  ::          liq_frac


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  Internal variables
!  ------------------

real, dimension(size(tau,1),size(tau,2),size(tau,3))   :: wgt_ice, wgt_liq
real, dimension(size(tau,1),size(tau,2),size(tau,3))   :: cwp, lwp, iwp

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!



    !---------------   COMPUTE OPTICAL DEPTH ---------------------------!



        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately


!       by Tony Gordon's v197 scheme

        WHERE (tempcld(:,:,:) .le. tk_all_ice)
               wgt_liq(:,:,:) = 0.
               wgt_ice(:,:,:) = 1.
        END WHERE

        WHERE (tempcld(:,:,:) .ge. tk_all_liq)
               wgt_liq(:,:,:) = 1.
               wgt_ice(:,:,:) = 0.
        END WHERE

        WHERE (tempcld(:,:,:) .gt. tk_all_ice .and. tempcld(:,:,:) &
               .lt. tk_all_liq)
               wgt_liq(:,:,:) = (tempcld(:,:,:) - tk_all_ice) / &
                                (tk_all_liq - tk_all_ice)
               wgt_ice(:,:,:) = 1. - wgt_liq(:,:,:)
        END WHERE
                     
        CWP(:,:,:) = tau(:,:,:,1) / &
                     (k_sw_liq * wgt_liq(:,:,:) + &
                      k_sw_ice * wgt_ice(:,:,:) )

        LWP(:,:,:) = wgt_liq(:,:,:) * CWP(:,:,:)
        IWP(:,:,:) = wgt_ice(:,:,:) * CWP(:,:,:)
liq_frac = wgt_liq



end subroutine CLOUD_OPT_PROP_tg2


!#######################################################################

  SUBROUTINE DIAG_CLOUD_RAD_INIT(do_crad_init)

!=======================================================================
! ***** INITIALIZE Predicted Cloud Scheme
!=======================================================================

!---------------------------------------------------------------------
! Argument (Intent inout)
!  do_crad_init - logical switch to be set = .true. after init is done
!---------------------------------------------------------------------
 logical, intent(inout) :: do_crad_init

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------
 integer             :: unit, io, logunit, ierr

!=====================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=diag_cloud_rad_nml, iostat=io)
  ierr = check_nml_error(io,"diag_cloud_rad_nml")
#else
  if( FILE_EXIST( 'input.nml' ) ) then
! -------------------------------------
         unit = open_namelist_file ()
   io = 1
   do while( io .ne. 0 )
      READ ( unit,  nml = diag_cloud_rad_nml, iostat = io, end = 10 ) 
      ierr = check_nml_error(io,'diag_cloud_rad_nml')
   end do
10 continue
   call close_file (unit)
! -------------------------------------
  end if
#endif

!   **** call cloud_rad_init to read namelist containing L2STREM  ****
       call cloud_rad_init()

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit, nml=diag_cloud_rad_nml)
      endif

!-------------------------------------------------------------------
  do_crad_init = .true.
  module_is_initialized = .true.
END SUBROUTINE DIAG_CLOUD_RAD_INIT

SUBROUTINE DIAG_CLOUD_RAD_END

  module_is_initialized = .false.

END SUBROUTINE DIAG_CLOUD_RAD_END


end MODULE DIAG_CLOUD_RAD_MOD



