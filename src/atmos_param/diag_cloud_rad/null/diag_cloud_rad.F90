MODULE DIAG_CLOUD_RAD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!CLOUD RADIATIVE PROPERTIES
!
!       May-Oct  1998 -> Sep 2000
!       Contact persons: Tony Gordon, Bill Stern (for modified code)
!                        Steve Klein (for original Fotran 90 code)
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

use       fms_mod, only:  error_mesg, FATAL, file_exist,    &
                          check_nml_error, open_namelist_file,       &
                          mpp_pe, mpp_root_pe, close_file, &
                          write_version_number, stdlog

! Steve Klein's Cloud_Rad module
  use Cloud_Rad_Mod, ONLY: CLOUD_RAD, CLOUD_RAD_INIT,              &
                          cloud_rad_k_diag


!-------------------------------------------------------------------
 
  implicit none

!-------------------------------------------------------------------

!        The module contains the following:
!
!SUBROUTINES
!
!            CLOUD_TAU_DRIVER
!                        calls a sequence of suboutines to compute
! cloud optical and radiative properties, as detailed
!                        below -->
!            CLOUD_PRES_THICK_FOR_TAU
!                        computes cloud-type dependent set of pressure 
!                        thicknesses for each distinct cloud layer, which are  
!                        used to parameterize cloud optical depths
!            CLOUD_OPTICAL_DEPTHS
! Specify / crudely parameterize cloud optical depths
! for distinct cloud layers,incorporating a 
! parameterization scheme for non-anvil cirrus
!                        proposed by Harshvardhan, based upon observations by
!                        Platt and Harshvardhan
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!          

 private



!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: diag_cloud_rad.F90,v 17.0 2009/07/21 02:54:14 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------


!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 public cloud_tau_driver, diag_cloud_rad_init, diag_cloud_rad_end, &
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       ------
!INPUT:
!       ------
!
!      TEMPCLD    cloud layer mean temperature (degrees Kelvin) of distinct
!                    cloud layers
!      COSZEN     cosine of the zenith angle


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

real, intent(in), dimension(:,:,:)   :: tempcld
real, intent(in), dimension(:,:)     :: qmix_kx
real, intent(in), dimension(:,:)     :: coszen
real, intent(in), dimension(:,:,:,:) :: tau


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real, intent (out), dimension(:,:,:), optional  :: r_uv,r_nir,  &
                                                   ab_uv,ab_nir,em_lw

!  *****************************************************************



!--------------------------------------------------------------------

      call error_mesg('CLOUD_TAU_DRIVER', &
      'This module is not supported as part of the public release', FATAL)

END SUBROUTINE CLOUD_TAU_DRIVER





!########################################################################
subroutine cloud_pres_thick_for_tau (nclds,icld,cldtop,cldbas, &
     &          delp_true,lhight,lhighb, lmidt, lmidb, llowt,lk,delp, &
     &          phalf, psfc )

! This subroutine calculates a special cloud-type dependent set of 
! cloud pressure thicknesses that will be employed to compute cloud optical
! depths.

!===================================================================

! Arguments (intent in)

real,    intent(in), dimension(:,:,:) :: phalf, delp_true
integer, intent(in), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(in), dimension(:,:)   :: nclds
integer, intent(in), dimension(:,:)   :: lhight, lhighb, lmidt, lmidb, llowt
integer, intent(in)                   :: lk
real,    intent(in), dimension(:,:)   :: psfc

!      INPUT
!      ------

!       PHALF      pressure at half levels (Pascals)
!                    NOTE: it is assumed that phalf(j+1) > phalf(j)
!       PSFC       Surface pressure field
!       NCLDS      number of (random overlapping) clouds in column and also
!                    the current # for clouds to be operating on
!       ICLD       marker array of cloud types/heights (at cloud levels)
!       CLDTOP     index of cloud tops (at cloud levels)
!       CLDBAS     index of cloud bottoms (at cloud levels)
!       DELP_TRUE  true cloud pressure thickness of distinct cloud layers 
!                    (at cloud levels - in Pa)
!       LHIGHT     vertical level index upper limit for high cloud tops
!       LHIGHB     vertical level index lower limit for high cloud bases
!       LMIDT      vertical level index upper limit for mid cloud tops
!       LMIDB      vertical level index lower limit for mid cloud bases
!       LLOWT      vertical level index upper limit for low cloud tops
!       LK         vertical level below which no low cloud bases can exist

!===================================================================

! Arguments (intent out)

real, intent(out), dimension(:,:,:) :: delp

!      OUTPUT
!      ------

!      DELP     cloud pressure thickness used to calculate cloud optical depths 
!                   of distinct cloud layers
!=======================================================================




      call error_mesg('cloud_pres_thick_for_tau', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cloud_pres_thick_for_tau

!########################################################################

subroutine cloud_optical_depths2 (nclds,icld,cldtop,cldbas,tempcld,delp, &
     &          tau,phalf,           liq_frac )

! This subroutine specifies/crudely parameterizes cloud optical depths 
! of non-anvil cirrus clouds based upon a parameterization scheme 
! resembling that proposed by Harshvardhan,
! (which is based on observations by Platt and Harshvardhan).

!===================================================================

! Arguments (intent in)

real,    intent(in), dimension(:,:,:) :: phalf, delp, tempcld
integer, intent(in), dimension(:,:,:) :: cldtop, cldbas, icld
integer, intent(in), dimension(:,:)   :: nclds

!      INPUT
!      ------

!      PHALF      pressure at half levels (Pascals)
!                   NOTE: it is assumed that phalf(j+1) > phalf(j)
!      NCLDS      number of (random overlapping) clouds in column and also
!                   the current # for clouds to be operating on
!      ICLD       marker array of cloud types (at cloud levels)
!      CLDTOP     index of cloud tops (at cloud levels)
!      CLDBAS     index of cloud bottoms (at cloud levels)
!      TEMPCLD    cloud layer mean temperature 
!                   (degrees Kelvin, at cloud levels)
!      DELP       cloud pressure thickness used for cloud optical depth 
!                   (at cloud levels)

!===================================================================

! Arguments (intent out)

real, intent(inout), dimension(:,:,:,:) :: tau
real, intent(inout), dimension(:,:,:)   :: liq_frac

!      OUTPUT
!      ------

!      TAU     cloud optical depth (at cloud levels)

!=======================================================================



      call error_mesg('cloud_optical_depths2', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cloud_optical_depths2

subroutine cloud_optical_depths (nclds,icld,cldtop,cldbas,tempcld,delp, &
     &          tau,phalf )

! This subroutine specifies/crudely parameterizes cloud optical depths 
! of non-anvil cirrus clouds based upon a parameterization scheme 
! resembling that proposed by Harshvardhan,
! (which is based on observations by Platt and Harshvardhan).

!===================================================================

! Arguments (intent in)

real,    intent(in), dimension(:,:,:) :: phalf, delp, tempcld
integer, intent(in), dimension(:,:,:) :: cldtop, cldbas ,icld
integer, intent(in), dimension(:,:)   :: nclds

!      INPUT
!      ------

!      PHALF      pressure at half levels (Pascals)
!                   NOTE: it is assumed that phalf(j+1) > phalf(j)
!      NCLDS      number of (random overlapping) clouds in column and also
!                   the current # for clouds to be operating on
!      ICLD       marker array of cloud types (at cloud levels)
!      CLDTOP     index of cloud tops (at cloud levels)
!      CLDBAS     index of cloud bottoms (at cloud levels)
!      TEMPCLD    cloud layer mean temperature 
!                   (degrees Kelvin, at cloud levels)
!      DELP       cloud pressure thickness used for cloud optical depth 
!                   (at cloud levels)

!===================================================================

! Arguments (intent out)

real, intent(out), dimension(:,:,:,:) :: tau

!      OUTPUT
!      ------

!      TAU     cloud optical depth (at cloud levels)

!=======================================================================

      call error_mesg('cloud_optical_depths', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cloud_optical_depths

!########################################################################

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real, intent (in),  dimension(:,:,:)   :: liq_frac
real, intent (in),  dimension(:,:,:,:) :: tau                 
real, intent (out), dimension(:,:,:)   :: em_lw

      call error_mesg('CLOUD_OPT_PROP_tg_lw', &
      'This module is not supported as part of the public release', FATAL)

end subroutine CLOUD_OPT_PROP_tg_lw

subroutine CLOUD_OPT_PROP_tg_sw( liq_frac,       &
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
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!
real,     intent (in),  dimension(:,:)     :: qmix_kx, cosz
logical,  intent (in),  dimension(:,:,:)   :: direct           
real,     intent (in),  dimension(:,:,:,:) :: tau
real,     intent (in),  dimension(:,:,:)   :: liq_frac
real,     intent (out), dimension(:,:,:)   :: cuvrf, cirrf,   &
                                              cuvab, cirab

      call error_mesg('CLOUD_OPT_PROP_tg_sw', &
      'This module is not supported as part of the public release', FATAL)

end subroutine CLOUD_OPT_PROP_tg_sw




subroutine CLOUD_OPT_PROP_tg2 (tau, tempcld, liq_frac)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      tempcld      cloud layer mean temperature (degrees Kelvin), with
!                   compressed cloud layer index.
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real, intent (in),  dimension(:,:,:,:) :: tau
real, intent (in),  dimension(:,:,:)   :: tempcld

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

real, intent (out), dimension(:,:,:)   :: liq_frac


      call error_mesg('CLOUD_OPT_PROP_tg2', &
      'This module is not supported as part of the public release', FATAL)

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
! --- Output version
!---------------------------------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

!-------------------------------------------------------------------
      module_is_initialized = .true.

      call error_mesg('DIAG_CLOUD_RAD_INIT', &
      'This module is not supported as part of the public release', FATAL)

END SUBROUTINE DIAG_CLOUD_RAD_INIT

!#######################################################################

SUBROUTINE DIAG_CLOUD_RAD_END

!-------------------------------------------------------------------

  module_is_initialized = .false.


      call error_mesg('DIAG_CLOUD_RAD_END', &
      'This module is not supported as part of the public release', FATAL)

END SUBROUTINE DIAG_CLOUD_RAD_END

end MODULE DIAG_CLOUD_RAD_MOD



