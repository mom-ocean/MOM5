
                        MODULE LONGWAVE_MOD

!-----------------------------------------------------------------------

      USE RDPARM_MOD, ONLY: LMAX
      USE RDPARM_MOD, ONLY: LM1,LP1,LP2,LL,LLP1,LLM1,LP1M,LP1V,LL3P
      USE RDPARM_MOD, ONLY: NBLW,NBLX,NBLY,NBLM,INLTE,INLTEP,NNLTE

      USE HCONST_MOD, ONLY: DIFFCTR,GINV,P0,P0INV,GP0INV,P0XZP2,P0XZP8
      USE HCONST_MOD, ONLY: RADCON,RADCON1,RATH2OMW,SECPDA

      Use    FMS_Mod, ONLY:  Error_Mesg, FATAL, NOTE, mpp_pe, &
                             mpp_root_pe, write_version_number

      Use CO2_Data_Mod, ONLY:  CO251,CO258,CDT51,CDT58,C2D51,C2D58, &
                               CO2M51,CO2M58,CDTM51,CDTM58,C2DM51,  &
                               C2DM58, STEMP,GTEMP, B0,B1,B2,B3
      Use CO2_Data_Mod, ONLY:  CO231,CO238,CDT31,CDT38,C2D31,C2D38
      Use CO2_Data_Mod, ONLY:  CO271,CO278,CDT71,CDT78,C2D71,C2D78
      Use CO2_Data_Mod, ONLY:  CO211,CO218


!     -----------------------------------------------------------
implicit none
private

!-----------------------------------------------------------------------
!--------------------- G L O B A L   D A T A ---------------------------
!-----------------------------------------------------------------------
!
!    Random band parameters for the longwave calcualtions using
!    10 cm-1 wide bands. The 15 um co2 complex is 2 bands,
!    560-670 and 670-800 cm-1. Ozone coefficients are in 3 bands,
!    670-800 (14.1 um), 990-1070 and 1070-1200 (9.6 um).
!    The (NBLW) bands now include: 
!
!                56 BANDS, 10  CM-1 WIDE    0  -   560  CM-1
!                 2 BANDS, 15 UM COMPLEX  560  -   670  CM-1
!                                         670  -   800  CM-1
!                 3 "CONTINUUM" BANDS     800  -   900  CM-1
!                                         900  -   990  CM-1
!                                        1070  -   1200 CM-1
!                 1 BAND FOR 9.6 UM BAND  990  -   1070 CM-1
!               100 BANDS, 10 CM-1 WIDE  1200  -   2200 CM-1
!                 1 BAND FOR 4.3 UM SRC  2270  -   2380 CM-1
!
!    Thus NBLW presently equals    163
!    All bands are arranged in order of increasing wavenumbers.
! 
!      ARNDM   =   RANDOM "A" PARAMETER FOR (NBLW) BANDS
!      BRNDM   =   RANDOM "B" PARAMETER FOR (NBLW) BANDS
!      BETAD   =   CONTINUUM COEFFICIENTS FOR (NBLW) BANDS
!      AP,BP   =   CAPPHI COEFFICIENTS FOR (NBLW) BANDS 
!      ATP,BTP =   CAPPSI COEFFICIENTS FOR (NBLW) BANDS 
!      BANDLO  =   LOWEST FREQUENCY IN EACH OF (NBLW) FREQ. BANDS 
!      BANDHI  =   HIGHEST FREQUENCY IN EACH OF (NBLW) FREQ. BANDS
!      AO3RND  =   RANDOM "A" PARAMETER FOR OZONE IN (3) OZONE
!                  BANDS.
!      BO3RND  =   RANDOM "B" PARAMETER FOR OZONE IN (3) OZONE
!                  BANDS
!      AB15    =   THE PRODUCT ARNDM*BRNDM FOR THE TWO BANDS
!                  REPRESENTING THE 15 UM BAND COMPLEX OF CO2 
!
!     Data for ARNDM,BRNDM,AP,BP,ATP,BTP,AO3RND,BO3RND are obtained
!     by using the AFGL 1982 catalog. Continuum coefficients are from
!     Roberts (1976). This data was formerly in COMMON /BANDTA/.

      REAL  ARNDM(NBLW),BRNDM(NBLW),BETAD(NBLW),AP(NBLW), &
            BP(NBLW),ATP(NBLW),BTP(NBLW),BANDLO(NBLW),    &
            BANDHI(NBLW),AO3RND(3),BO3RND(3),AB15(2)

!-----------------------------------------------------------------------
!
!    Random band parameters for the longwave calculations using
!    comboned wide frequency bands between 160 and 1200 cm-1, 
!    as well as the 2270-2380 band for source calculations.
!
!        BANDS 1-8: COMBINED WIDE FREQUENCY BANDS FOR 160-560 CM-1
!        BANDS 9-14: FREQUENCY BANDS,AS IN BANDTA (NARROW BANDS)
!                    FOR 560-1200 CM-1
!        BAND  15:  FREQUENCY BAND 2270-2380 CM-1,USED FOR SOURCE 
!                   CALCULATION ONLY
!
!        Thus NBLY presently equals   15
! 
!        Bands are arranged in order of increasing wavenumber.
!
!      ACOMB       =   RANDOM "A" PARAMETER FOR (NBLY) BANDS
!      BCOMB       =   RANDOM "B" PARAMETER FOR (NBLY) BANDS
!      BETACM      =   CONTINUUM COEFFICIENTS FOR (NBLY) BANDS
!      APCM,BPCM   =   CAPPHI COEFFICIENTS FOR (NBLY) BANDS 
!      ATPCM,BTPCM =   CAPPSI COEFFICIENTS FOR (NBLY) BANDS 
!      BDLOCM      =   LOWEST FREQUENCY IN EACH OF (NBLY) FREQ. BANDS 
!      BDHICM      =   HIGHEST FREQUENCY IN EACH OF (NBLY) FREQ. BANDS
!      AO3CM       =   RANDOM "A" PARAMETER FOR OZONE IN (3) OZONE
!                      BANDS.
!      BO3CM       =   RANDOM "B" PARAMETER FOR OZONE IN (3) OZONE
!                      BANDS
!      AB15CM      =   THE PRODUCT ARNDM*BRNDM FOR THE TWO BANDS
!                      REPRESENTING THE 15 UM BAND COMPLEX OF CO2 
!      BETINC      =   CONT.COEFFICIENT FOR A SPECIFIED WIDE
!                      FREQ.BAND (800-990 AND 1070-1200 CM-1).
!      IBAND       =   INDEX NO OF THE 40 WIDE BANDS USED IN
!                      COMBINED WIDE BAND CALCULATIONS. IN OTHER
!                      WORDS,INDEX TELLING WHICH OF THE 40 WIDE 
!                      BANDS BETWEEN 160-560 CM-1 ARE INCLUDED IN 
!                      EACH OF THE FIRST 8 COMBINED WIDE BANDS
!
!     Data for ACOMB,BCOMB,APCM,BPCM,ATPCM,BTPCM,AO3CM,BO3CM are
!     obtained by using the AFGL 1982 catalog. Continuum coefficients 
!     are from Roberts (1976). IBAND index values are obtained by 
!     experimentation. This data was formerly in COMMON /BDCOMB/.

      INTEGER  IBAND(40)
      REAL  ACOMB(NBLY),BCOMB(NBLY),                        &
            BETACM(NBLY),APCM(NBLY),BPCM(NBLY),ATPCM(NBLY), &
            BTPCM(NBLY),BDLOCM(NBLY),BDHICM(NBLY),BETINC,   &
            AO3CM(3),BO3CM(3),AB15CM(2) 


!-----------------------------------------------------------------------
!
!    Random band parameters for specific wide bands. At present,
!    the information consists of:  1) random model parameters for
!    the 15 um band, 560-800 cm-1; 2) the continuum coefficient for
!    the 800-990, 1070-1200 cm-1 band.
!
!    specifically:  
!      AWIDE       =   RANDOM "A" PARAMETER FOR  BAND 
!      BWIDE       =   RANDOM "B" PARAMETER FOR  BAND 
!      BETAWD      =   CONTINUUM COEFFICIENTS FOR BAND
!      APWD,BPWD   =   CAPPHI COEFFICIENTS FOR  BAND
!      ATPWD,BTPWD =   CAPPSI COEFFICIENTS FOR BAND 
!      BDLOWD      =   LOWEST FREQUENCY IN EACH  FREQ  BAND 
!      BDHIWD      =   HIGHEST FREQUENCY IN EACH FREQ  BAND 
!      AB15WD      =   THE PRODUCT ARNDM*BRNDM FOR THE ONE BAND 
!                      REPRESENTING THE 15 UM BAND COMPLEX OF CO2 
!      BETINW      =   CONT.COEFFICIENT FOR A SPECIFIED WIDE
!                      FREQ.BAND (800-990 AND 1070-1200 CM-1).
!      SKO2D       =   1./BETINW, USED IN SPA88 FOR CONT. COEFFS
!      SKC1R       =   BETAWD/BETINW, USED FOR CONT. COEFF. FOR 
!                      15 UM BAND IN FST88
!      SKO3R       =   RATIO OF CONT. COEFF. FOR 9.9 UM BAND TO 
!                        BETINW, USED FOR 9.6 UM CONT COEFF IN FST88
!
!     Data for AWIDE,BWIDE,APWD,BPWD,ATPWD,BTPWD,AO3WD,BO3WD are
!     obtained by using the AFGL 1982 catalog. Continuum coefficients 
!     are from Roberts (1976). This data was formerly in
!     COMMON /BDWIDE/.

      REAL  AWIDE,BWIDE,BETAWD,    &
            APWD,BPWD,ATPWD,BTPWD, &
            BDLOWD,BDHIWD,BETINW,  &
            AB15WD,SKO2D,SKC1R,SKO3R


!-----------------------------------------------------------------------
!
!       CLDFAC     =  CLOUD TRANSMISSION FUNCTION,ASSUMING RANDOM
!                       OVERLAP

      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CLDFAC

!-----------------------------------------------------------------------
!
!     Basic quantities computed in SUBROUTINE LWRAD and used in
!     the remaining longwave routines (formally COMMON /KDACOM/):
!
!       QH2O     =  H2O MASS MIXING RATIO,MULTIPLIED BY THE 
!                     DIFFUSIVITY FACTOR (DIFFCTR)
!       P        =  PRESSURE AT FLUX LEVELS OF MODEL
!       DELP2    =  PRESSURE DIFFERENCE BETWEEN FLUX LEVELS 
!       DELP     =  INVERSE OF DELP2
!       TTTT     =  TEMPERATURE ASSIGNED TO MODEL FLUX LEVELS 
!       VAR1     =  H2O OPTICAL PATH IN MODEL LAYERS (BETWEEN 
!                     FLUX LEVELS)
!       VAR2     =  PRESSURE-WEIGHTED H2O OPTICAL PATH IN MODEL LAYERS
!       VAR3     =  O3 OPTICAL PATH IN MODEL LAYERS 
!       VAR4     =  PRESSURE-WEIGHTED O3 OPTICAL PATH IN MODEL LAYERS 
!       CNTVAL   =  H2O CONTINUUM PATH IN MODEL LAYERS FOR THE
!                     800-990 AND 1070-1200 CM-1 COMBINED BAND

      REAL, ALLOCATABLE, DIMENSION(:,:) :: QH2O,P,DELP2,DELP,TTTT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VAR1,VAR2,VAR3,VAR4,CNTVAL

!-----------------------------------------------------------------------
!
!     Flux quantities computed by the radiation code, used for
!     diagnostic purposes (formally COMMON /RDFLUX/): 
!
!       FLX1E1     =  FLUX AT TOP FOR 0-160,1200-2200 CM-1 RANGE
!       GXCTS      =  FLUX AT TOP FOR 160-1200 CM-1 RANGE 
!       FCTSG      =  CTS FLUX AT GROUND. USED TO OBTAIN GXCTS
!                              BY BANDS.

      REAL, ALLOCATABLE, DIMENSION(:)   :: FLX1E1,GXCTS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FCTSG

!-----------------------------------------------------------------------
!
!     Planck function values used for the radiative calculations
!     (formally COMMON /SRCCOM/):
!
!       SORC     =  PLANCK FCTN, AT MODEL TEMPERATURES, FOR ALL BANDS 
!                     USED IN CTS CALCULATIONS
!       CSOUR1   =  PLANCK FCTN FOR 560-670 CM-1 BAND 
!       CSOUR2   =  PLANCK FCTN FOR 670-800 CM-1 BAND 
!       CSOUR    =  PLANCK FCTN FOR 560-800 CM-1 BANDS
!       OSOUR    =  PLANCK FCTN FOR 990-1070 CM-1 BAND
!       SS1      =  PLANCK FCTN FOR 800-990,1070-1200 CM-1 BANDS

      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SORC
      REAL, ALLOCATABLE, DIMENSION(:,:) :: CSOUR1,CSOUR2,OSOUR,CSOUR,SS1

!-----------------------------------------------------------------------
!
!     Quantities precomputed in subroutine TABLE for use in
!     the longwave radiation module (formally COMMON /TABCOM/):
!
!        EM1     =  E1 FUNCTION, EVALUATED OVER THE 0-560 AND 
!                   1200-2200 CM-1 INTERVALS
!        EM1WDE  =  E1 FUNCTION, EVALUATED OVER THE 160-560 CM-1
!                   INTERVAL
!        TABLE1  =  E2 FUNCTION, EVALUATED OVER THE 0-560 AND 
!                   1200-2200 CM-1 INTERVALS
!        TABLE2  =  TEMPERATURE DERIVATIVE OF TABLE1
!        TABLE3  =  MASS DERIVATIVE OF TABLE1 
!        EM3     =  E3 FUNCTION, EVALUATED OVER THE 0-560 AND 
!                   1200-2200 CM-1 INTERVALS
!        SOURCE  =  PLANCK FUNCTION, EVALUATED AT SPECIFIED TEMPS. FOR
!                   BANDS USED IN CTS CALCULATIONS
!        DSRCE   =  TEMPERATURE DERIVATIVE OF SOURCE
!        INDX2   =  INDEX VALUES USED IN OBTAINING "LOWER TRIANGLE" 
!                   ELEMENTS OF AVEPHI,ETC.,IN FST88
!        KMAXV   =  INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" 
!                   ELEMENTS OF AVEPHI,ETC.,IN FST88
!        KMAXVM  =  KMAXV(LMAX),USED FOR DO LOOP INDICES 

      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDX1,INDX2,KMAXV

      INTEGER :: KMAXVM

      REAL, DIMENSION(28,NBLY) :: SOURCE,DSRCE

!     REAL  EM1   (28,180),EM1WDE(28,180),TABLE1(28,180),
!     COMMON /TABCOM/  EM1   (28,180),EM1WDE(28,180),TABLE1(28,180),
!    &                 TABLE2(28,180),TABLE3(28,180),EM3   (28,180)
      
!-----------------------------------------------------------------------
!
!     Transmission functions used for radiative computations, and
!     output heating rates and fluxes, except those needed out of
!     the radiative module (formally COMMON /TFCOM/): 
!
!       TO3      =  TRANSMISSION FCTN FOR THE 990-1070 CM-1 BAND
!                     O3(9.6 UM) + H2O CONTINUUM (NO LINES) 
!       CO21     =  TRANSMISSION FCTN FOR THE 560-800 CM-1 BAND 
!                     (AS 1 BAND). INCLUDES CO2 (IN LWRAD) AND
!                      H2O(L+C) AFTER MULTIPLICATION WITH "OVER"
!                      IN FST88 
!       EMISS    =  E2 EMISSIVITYY FCTN FOR H2O LINES (0-560,1200-2200
!                      CM-1). OBTAINED IN E1E288. 
!       EMISS2   =  TRANSMISSION FCTN FOR H2O CONTINUUM IN THE 800-990
!                      AND 1070-1200 CM-1 REGION, TAKEN AS 1 BAND 
!       AVEPHI   =  H2O OPTICAL PATHS BET. FLUX PRESSURES: INPUT TO 
!                      EMISSIVITY CALCULATIONS. 
!       TTEMP    =  TEMPERATURES USED AS INPUT FOR EMISSIVITY CALCS.
!       CTS      =  APPROX CTS HEATING RATES FOR 160-560 AND 800-990, 
!                      1070-1200 CM-1 RANGES
!       CTSO3    =  APPROX CTS HEATING RATES FOR 560-800,990-1070 CM-1
!                      RANGES 
!       EXCTS    =  EXACT CTS HEATING RATES FOR 160-1200 CM-1 RANGE 
!       EXCTSN   =  EXACT CTS HEATING RATES, BY BANDS 
!       E1FLX    =  E1 EMISSIVITY FCTN FOR H2O LINES (0-560,1200-CM-1)
!       CO2NBL   =  CO2 TRANS. FCTNS. (NOT PRESSURE-INTEGRATED) FOR 
!                      ADJACENT LEVELS,OVER THE 560-800 CM-1 RANGE. 
!       CO2SP1   =  CO2 TRANS. FCTNS. (NOT PRESSURE-INTEGRATED) BET.
!                      A FLUX LEVEL AND SPACE, FOR THE 560-670 CM-1 
!                      RANGE. USED FOR EXACT CTS CALCS. 
!       CO2SP2   =  SAME AS CO2SP1, BUT FOR THE 670-800 CM-1 RANGE. 
!       CO2SP    =  SAME AS CO2SP1, BUT FOR THE 560-800 CM-1 BAND.
!                      USED FOR APPROX CTS CALCS. 
!       TO3SPC   =  O3 OPTICAL DEPTHS BET. A LEVEL AND SPACE. USED FOR
!                      EXACT CTS CALCS. 
!       TOTVO2   =  H2O CONTINUUM OPTICAL PATHS BET. SPACE AND A
!                      LEVEL, USING THE CNT. COEFFICIENT FOR THE
!                      1-BAND 800-990,1070-1200 CM-1 BAND. USED FOR 
!                      CTS CALCS. 

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: TO3,CO21,EMISS,EMISS2,AVEPHI
      REAL,ALLOCATABLE,DIMENSION(:,:)   :: CTS,CTSO3,EXCTS
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: EXCTSN
      REAL,ALLOCATABLE,DIMENSION(:,:)   :: E1FLX,CO2NBL,CO2SP1,CO2SP2
      REAL,ALLOCATABLE,DIMENSION(:,:)   :: CO2SP,TO3SPC,TOTVO2

!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: longwave.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized = .false.

!-----------------------------------------------------------------------

public LWRad, Rad_DeAlloc, longwave_init, longwave_end
public OSOUR, CSOUR, SS1, FLX1E1, GXCTS, FCTSG, CLDFAC, DELP2, DELP, &
       TO3, CO21, EMISS, EMISS2, CTS, EXCTS, EXCTSN, E1FLX, CO2SP,   &
       IBAND, BANDLO, BANDHI

      CONTAINS

!#######################################################################
!#######################################################################
      Subroutine longwave_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('longwave_init', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine longwave_init

!#######################################################################
!#######################################################################

      Subroutine longwave_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('longwave_end', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine longwave_end

!#######################################################################
!#######################################################################

      SUBROUTINE LWRAD (KTOP,KBTM,NCLDS,EMCLD,PRES,TEMP,RH2O,QO3,CAMT, &
                        RRVCO2,  HEATRA,GRNFLX,TOPFLX, LSFC,PSFC)

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)   :: NCLDS
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRES,TEMP,RH2O,QO3,CAMT
      REAL,    INTENT(IN)                   :: RRVCO2

      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: HEATRA
      REAL,   INTENT(OUT), DIMENSION(:,:)   :: GRNFLX,TOPFLX

      INTEGER, INTENT(IN),OPTIONAL, DIMENSION(:,:)   :: LSFC
      REAL   , INTENT(IN),OPTIONAL, DIMENSION(:,:)   :: PSFC

!---------------------------------------------------------------------

      call error_mesg('LWRAD', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE LWRAD

!#######################################################################

      SUBROUTINE RAD_DEALLOC

!---------------------------------------------------------------------

      call error_mesg('RAD_DEALLOC', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE RAD_DEALLOC

!#######################################################################

                     END MODULE LONGWAVE_MOD

