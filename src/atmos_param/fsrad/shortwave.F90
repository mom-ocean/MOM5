
                        MODULE SHORTWAVE_MOD

!-----------------------------------------------------------------------

      USE  RDPARM_MOD, ONLY: LMAX,LP1,LLP1,LP2,LLP2,NB
      USE  HCONST_MOD, ONLY: DIFFCTR,GINV,O3DIFCTR,RADCON

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

!implicit none
private

!------- interfaces -------
      PUBLIC  SWRAD, SHORTWAVE_INIT, SHORTWAVE_END

      character(len=128) :: version = '$Id: shortwave.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical            :: module_is_initialized = .false.

      integer :: IMAX
      CONTAINS

!#######################################################################

      SUBROUTINE SWRAD (NCLDS,KTOPSW,KBTMSW,PRESS,RH2O,QO3,CAMT, &
                        CUVRF,CIRRF,CIRAB,RRCO2,COSZRO,SSOLAR, &
                        SALB, FSW,DFSW,UFSW,HSW, LSFC,PSFC)

!-----------------------------------------------------------------------
!              WRAPPER FOR  SHORT WAVE RADIATION CODE
!     inserts surface albedo into appropriate cloud property arrays
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN), DIMENSION(:,:)   :: NCLDS
      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRESS,RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: CAMT,CUVRF,CIRRF,CIRAB
      REAL,    INTENT(IN)                   :: RRCO2
      REAL,    INTENT(IN), DIMENSION(:,:)   :: COSZRO,SSOLAR
      REAL,    INTENT(IN), DIMENSION(:,:)   :: SALB

      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   :: LSFC
         REAL, INTENT(IN), OPTIONAL, DIMENSION(:,:)   :: PSFC

!-----------------------------------------------------------------------
      REAL,DIMENSION(SIZE(CUVRF,1),SIZE(CUVRF,3)) :: CUVRF2,CIRRF2
      REAL,DIMENSION(SIZE(CAMT ,1),SIZE(CAMT ,3)) :: CAMT2

      INTEGER  i,j
!-----------------------------------------------------------------------

      IMAX=SIZE(PRESS,1)

      DO j=1,SIZE(PRESS,2)

      CUVRF2(:,:)=CUVRF(:,j,:)
      CIRRF2(:,:)=CIRRF(:,j,:)
      CAMT2 (:,:)=CAMT (:,j,:)

      DO i=1,IMAX
         IF (COSZRO(i,j) > 0.0) THEN
            CUVRF2(i,NCLDS(i,j)+2)=SALB(i,j)
            CIRRF2(i,NCLDS(i,j)+2)=SALB(i,j)
            CAMT2 (i,NCLDS(i,j)+2)=1.0
         ENDIF
      ENDDO

!----- check usage of optional arguments -------
      IOPT=0
      IF (PRESENT(LSFC)) IOPT=IOPT+1
      IF (PRESENT(PSFC)) IOPT=IOPT+2

!     ------------------
      SELECT CASE (IOPT)
!     ------------------
          CASE (0)
!     ------------------
      CALL SWRAD_ORIG (NCLDS(:,j),KTOPSW(:,j,:),KBTMSW(:,j,:), &
                       PRESS(:,j,:),RH2O(:,j,:),QO3(:,j,:),CAMT2(:,:), &
                       CUVRF2,CIRRF2,CIRAB(:,j,:),RRCO2, &
                       COSZRO(:,j),SSOLAR(:,j), &
                       FSW(:,j,:),DFSW(:,j,:),UFSW(:,j,:),HSW(:,j,:))
!     ------------------
          CASE (3)
!     ------------------
      CALL SWRAD_ORIG (NCLDS(:,j),KTOPSW(:,j,:),KBTMSW(:,j,:), &
                       PRESS(:,j,:),RH2O(:,j,:),QO3(:,j,:),CAMT2(:,:), &
                       CUVRF2,CIRRF2,CIRAB(:,j,:),RRCO2, &
                       COSZRO(:,j),SSOLAR(:,j), &
                       FSW(:,j,:),DFSW(:,j,:),UFSW(:,j,:),HSW(:,j,:), &
                       LSFC(:,j),PSFC(:,j))
!     ------------------
          CASE DEFAULT
!     ------------------
             Call Error_Mesg ('SWRAD in SHORTWAVE_MOD', &
                          'LSFC and PSFC must be used together.', FATAL)
!     ------------------
      END SELECT
!     ------------------

      ENDDO
!-----------------------------------------------------------------------

      END SUBROUTINE SWRAD
      
!#######################################################################

      SUBROUTINE SWRAD_ORIG (NCLDS,KTOPSW,KBTMSW,PRESS,RH2O,QO3,CAMT, &
                        CUVRF,CIRRF,CIRAB,RRCO2,COSZRO,SSOLAR, &
                        FSW,DFSW,UFSW,HSW,LSFC,PSFC)

!***********************************************************************
!                    SHORT WAVE RADIATION CODE
!***********************************************************************
!-----------------------------------------------------------------------
!                        INPUT PARAMETERS
!                        ----------------
!
!      NCLDS   =  NO. CLOUDS AT EACH GRID PT. 
!      KTOPSW  =  INDEX OF (FLUX LEVEL) PRESSURE OF CLOUD TOP, USED 
!                    IN THE SHORTWAVE PROGRAM 
!      KBTMSW  =  INDEX OF (FLUX LEVEL) PRESSURE OF CLOUD BOTTOM, 
!                    USED IN THE SHORTWAVE PROGRAM
!      PRESS   =  PRESSURE (CGS UNITS) AT DATA LEVELS OF MODEL
!      RH2O    =  MASS MIXING RATIO (G/G) OF H2O AT MODEL DATA LVLS.
!      QO3     =  MASS MIXING RATIO (G/G) OF O3 AT MODEL DATA LVLS. 
!      CAMT    =  CLOUD AMOUNTS OF CLOUDS (THEIR LOCATIONS ARE
!                    SPECIFIED IN THE KTOP/KBTM INDICES)
!      CUVRF   =  REFLECTIVITY OF CLOUDS IN THE VISIBLE FREQ. BAND
!                    USED IN SHORTWAVE CALCS. ONLY
!      CIRRF   =  REFLECTIVITY OF CLOUDS IN THE INFRARED FREQ. BAND 
!                    USED IN SHORTWAVE CALCS. ONLY
!      CIRAB   =  ABSORPTIVITY OF CLOUDS IN THE INFRARED FREQ. BAND 
!                    USED IN SHORTWAVE CALCS. ONLY
!      RRCO2   =  MASS MIXING RATIO (G/G) OF CO2,USED IN SHORTWAVE
!                    CALCS. ONLY (scalar)
!      COSZRO  =  ZENITH ANGLE AT GRID PT. USED ON SHORTWAVE CALCS. 
!      SSOLAR  =  TOTAL SOLAR FLUX EITHER FOR THE TIMESTEP OR AVERAGED
!                 OVER THE DAY OR YEAR,
!                 EQUALS THE SOLAR CONSTANT x NORMALIZED SOLAR FLUX
!                 (INCLUDES THE COS ZENITH ANGLE AND DAY FRACTION)
!                 (AT PRESENT,IN LY/MIN).
!
!      LSFC    =  Vertical index of the lowest model level,
!                    dimensioned by IMAX.
!      PSFC    =  Surface pressure
!
!-----------------------------------------------------------------------
!                        OUTPUT PARAMETERS
!                        -----------------
!
!      FSW     = NET RADIATION (UP-DOWN) IN CGS UNITS AT ALL
!                PRESSURE LEVELS
!     DFSW     = DOWNWARD RADIATION AT ALL PRESSURE LEVELS
!     UFSW     = UPWARD RADIATION AT ALL PRESSURE LEVELS
!      HSW     = SHORTWAVE HEATING RATES IN K/DAY FOR PRESSURE
!                LAYERS.
!
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN), DIMENSION(:)   :: NCLDS
      INTEGER, INTENT(IN), DIMENSION(:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:) :: PRESS,RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:) :: CAMT,CUVRF,CIRRF,CIRAB
      REAL,    INTENT(IN)                 :: RRCO2
      REAL,    INTENT(IN), DIMENSION(:)   :: COSZRO,SSOLAR

      REAL,   INTENT(OUT), DIMENSION(:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:)   :: LSFC
      REAL,    INTENT(IN), OPTIONAL, DIMENSION(:)   :: PSFC
       
!-----------------------------------------------------------------------
!                        D I M E N S I O N
!    &  NCLDS(IMAX),KTOPSW(IMAX,LP1),KBTMSW(IMAX,LP1)
!    2 ,PRESS(IMAX,LP1),RH2O(IMAX,LMAX),QO3(IMAX,LMAX)
!    3 ,CAMT(IMAX,LP1),CUVRF(IMAX,LP1),CIRRF(IMAX,LP1),CIRAB(IMAX,LP1)
!    4 ,COSZRO(IMAX),SSOLAR(IMAX),LSFC(IMAX),PSFC(IMAX)
!                        D I M E N S I O N
!    &  FSW(IMAX,LP1),DFSW(IMAX,LP1),UFSW(IMAX,LP1),HSW(IMAX,LMAX)
!-----------------------------------------------------------------------
!***********************************************************************
      LOGICAL BCLDS,BJTOP
!-----------------------------------------------------------------------
                         DIMENSION &
        BCLDS(IMAX,LP1),BJTOP(IMAX,LP1), &
        ICNT(LP1),ICNT1(LP1),IINCL(LP1),INDX4(IMAX,LP1),INDXK(IMAX), &
        IBETCL(LP1) 
!-----------------------------------------------------------------------
                         DIMENSION &
        DFN(IMAX,LP1,NB),UFN(IMAX,LP1,NB), &
        TTD(IMAX,LP1,NB),TTU(IMAX,LP1,NB), &
        PP    (IMAX,LP1),PPTOP (IMAX,LP1), &
        DP    (IMAX,LP1),DPCLD (IMAX,LP1), &
        CR    (IMAX,LP1),CT    (IMAX,LP1), &
        TDCL1 (IMAX,LP1),TDCL2 (IMAX,LP1),TUCL1 (IMAX,LP1), &
        TUCL1I(IMAX,LP1),TDCL1I(IMAX,LP1),TDCL2I(IMAX,LP1), &
        TCLU  (IMAX,LP1),TCLD  (IMAX,LP1),ALFAU (IMAX,LP1), &
        UFNCLU(IMAX,LP1),UFNCLD(IMAX,LP1), &
        DFNCLU(IMAX,LP1),DFNCLD(IMAX,LP1), &
        TEMP1(IMAX),TEMP2(IMAX),TEMP3(IMAX),TEMP4(IMAX), &
        TEMP5(IMAX),TEMP6(IMAX),ALFA (IMAX), &
        VV   (IMAX),REFL (IMAX),SECZ (IMAX),RRAY (IMAX)
!***********************************************************************
!-------DIMENSION OF VARIABLES EQUIVALENCED TO THOSE PREVIOUSLY---------
!***********************************************************************
!-----------------------------------------------------------------------
!                        D I M E N S I O N
      DIMENSION PPBOT(IMAX,LP1)
      DIMENSION DFNTOP(IMAX,NB)
      DIMENSION UD(IMAX,LP1),UR(IMAX,LP1)
!     DIMENSION UCO2(IMAX,LLP2),UO3(IMAX,LLP2),ACO2(IMAX,LLP2)
!     DIMENSION AO3(IMAX,LLP2)
      DIMENSION VTROW1(IMAX,LP1)
      DIMENSION VTROW2(IMAX,LP1),VTROW3(IMAX,LP1)
      DIMENSION FF(IMAX,LP1),FFCO2(IMAX,LP1),FFO3(IMAX,LP1)
      DIMENSION PR2(IMAX,LP1)
      DIMENSION DU(IMAX,LP1),DUCO2(IMAX,LP1),DUO3(IMAX,LP1)
      DIMENSION ADCO2(IMAX,LP1),AUCO2(IMAX,LP1),UDCO2(IMAX,LP1), URCO2(IMAX,LP1)
      DIMENSION ABSDO3(IMAX,LP1),ABSUO3(IMAX,LP1),UDO3(IMAX,LP1), URO3(IMAX,LP1)
      DIMENSION CR1D(IMAX*LP1),ALFU1D(IMAX*LP1),TCLU1D(IMAX*LP1), TCLD1D(IMAX*LP1)
!--DIMENSIONS OF LOCAL DATA VARIABLES---
      DIMENSION ABCFF(NB),PWTS(NB)
!***********************************************************************
!-----------------------------------------------------------------------
!     EQUIVALENCE (ADCO2(1,1),ACO2(1,1)),(AUCO2(1,1),ACO2(1,LP2))
!     EQUIVALENCE (UDCO2(1,1),UCO2(1,1)),(URCO2(1,1),UCO2(1,LP2))
!     EQUIVALENCE (ABSDO3(1,1),AO3(1,1)),(ABSUO3(1,1),AO3(1,LP2))
!     EQUIVALENCE (UDO3(1,1),UO3(1,1)),(URO3(1,1),UO3(1,LP2))

!     EQUIVALENCE (CR,UD),(CT,UR)
!     EQUIVALENCE (TDCL1I,ABSDO3,PPBOT),(TDCL2I,ABSUO3)
!     EQUIVALENCE (UFNCLU,UDO3),(UFNCLD,URO3)
!     EQUIVALENCE (DFNCLU,UDCO2),(DFNCLD,URCO2)
!     EQUIVALENCE (TCLU,ADCO2),(TCLD,AUCO2)
!     EQUIVALENCE (TTD(1,1,1),FF(1,1))
!     EQUIVALENCE (TTD(1,1,2),FFCO2(1,1))
!     EQUIVALENCE (TTD(1,1,3),PR2(1,1))
!     EQUIVALENCE (TTD(1,1,4),DU(1,1))
!     EQUIVALENCE (TTD(1,1,5),DUCO2(1,1))
!     EQUIVALENCE (TTD(1,1,6),DUO3(1,1))
!     EQUIVALENCE (TTD(1,1,7),VTROW1(1,1))
!     EQUIVALENCE (TTD(1,1,8),VTROW2(1,1))
!     EQUIVALENCE (TTD(1,1,9),VTROW3(1,1))
!     EQUIVALENCE (TTU(1,1,1),DFNTOP(1,1))
!     EQUIVALENCE (CR1D,CR),(ALFU1D,ALFAU),(TCLD1D,TCLD),(TCLU1D,TCLU)
!-----------------------------------------------------------------------
      DATA ABCFF /2*4.0E-5,.002,.035,.377,1.95,9.40,44.6,190./
      DATA PWTS /.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
      DATA CFCO2,CFO3 /508.96,466.64/
      DATA REFLO3 /1.9/
      DATA RRAYAV /0.144/
!-----------------------------------------------------------------------

      DO K=1,LP1
      DO I=1,IMAX
         FF   (I,K)=DIFFCTR
         FFCO2(I,K)=DIFFCTR
         FFO3 (I,K)=O3DIFCTR
      ENDDO
      ENDDO

!------ NOTE: converting pressures (PRESS) to CGS units -------

      DO K=2,LMAX
      DO I=1,IMAX
!CCC     PP(I,K)=0.50*(PRESS(I,K)+PRESS(I,K-1))
         PP(I,K)=5.0*(PRESS(I,K)+PRESS(I,K-1))
      ENDDO
      ENDDO

      DO I=1,IMAX
         SECZ(I)=35./SQRT(1224.*COSZRO(I)*COSZRO(I)+1.0)
!        SECZ(I)=1./COSZRO(I)
      ENDDO

      IF (.not.PRESENT(LSFC)) THEN
         DO I=1,IMAX
            PP(I,1)=0.00
            PP(I,LP1)=10.*PRESS(I,LP1)
         ENDDO
      ELSE
         DO I=1,IMAX
            PP(I,1)=0.00
            PP(I,LP1)=10.*PRESS(I,LP1)
            PP(I,LSFC(I)+1)=10.*PSFC(I)
         ENDDO
      ENDIF

      DO K=1,LMAX
      DO I=1,IMAX
         DP(I,K)=PP(I,K+1)-PP(I,K)
      ENDDO
      ENDDO

      IF (.not.PRESENT(LSFC)) THEN
         DO K=1,LMAX
         DO I=1,IMAX
!CCC        PR2(I,K)=0.50*(PP(I,K)+PP(I,K+1))/PRESS(I,LP1)
            PR2(I,K)=0.050*(PP(I,K)+PP(I,K+1))/PRESS(I,LP1)
         ENDDO
         ENDDO
      ELSE
         DO K=1,LMAX
         DO I=1,IMAX
!CCC        PR2(I,K)=0.50*(PP(I,K)+PP(I,K+1))/PSFC(I)
            PR2(I,K)=0.050*(PP(I,K)+PP(I,K+1))/PSFC(I)
         ENDDO
         ENDDO
      ENDIF

      DO K=1,LMAX
      DO I=1,IMAX
         DUO3 (I,K)=QO3 (I,K)*DP(I,K)*GINV
         DUCO2(I,K)=RRCO2    *DP(I,K)*GINV
         DU   (I,K)=RH2O(I,K)*DP(I,K)*GINV
      ENDDO
      ENDDO

      DO I=1,IMAX
         JTOP=KTOPSW(I,2) 
         DO K=1,JTOP 
            FFO3 (I,K)=SECZ(I) 
            FFCO2(I,K)=SECZ(I)
            FF   (I,K)=SECZ(I) 
         ENDDO
      ENDDO

!-----------------------------------------------------------------------
!     CALCULATE PRESSURE-WEIGHTED OPTICAL PATHS IN UNITS OF G/CM2.
!     PRESSURE WEIGHTING IS BY P**0.5 
!     UD IS THE DOWNWARD PATH,UR THE UPWARD PATH, 
!     AND THE CALCULATION IS MADE BY TAKING A PATH WITH AN ANGLE
!     OF (SECZ) FROM THE TOP OF THE ATMOSPHERE TO THE TOPMOST CLOUD,
!     THEN USING THE DIFFUSIVITY FACTOR (1.66) TO THE SURFACE AND FOR 
!     REFLECTED RADIATION. THE CODE BELOW REFLECTS THIS.
!-----------------------------------------------------------------------
!*****************************************
      IF (.not.PRESENT(LSFC)) THEN
         DO K=1,LMAX
         DO I=1,IMAX
            VTROW1(I,K)=DU   (I,K)*PR2(I,K)
            VTROW2(I,K)=DUCO2(I,K)*PR2(I,K)*CFCO2
            VTROW3(I,K)=DUO3 (I,K)*CFO3
         ENDDO
         ENDDO
      ELSE
         DO K=1,LMAX
         DO I=1,IMAX
            VTROW1(I,K)=0.00
            VTROW2(I,K)=0.00
            VTROW3(I,K)=0.00
         ENDDO
         ENDDO
         DO I=1,IMAX
            LMA=LSFC(I)
            DO K=1,LMA
               VTROW1(I,K)=DU   (I,K)*PR2(I,K)
               VTROW2(I,K)=DUCO2(I,K)*PR2(I,K)*CFCO2
               VTROW3(I,K)=DUO3 (I,K)*CFO3
            ENDDO
         ENDDO
      ENDIF
!*****************************************

      DO I=1,IMAX
         UD   (I,1)=0.00 
         UDCO2(I,1)=0.00
         UDO3 (I,1)=0.00 
      ENDDO

      DO K=2,LP1
!     DO I=1,IMAX
         UD   (:,K)=UD   (:,K-1)+VTROW1(:,K-1)*FF   (:,K) 
         UDCO2(:,K)=UDCO2(:,K-1)+VTROW2(:,K-1)*FFCO2(:,K)
         UDO3 (:,K)=UDO3 (:,K-1)+VTROW3(:,K-1)*FFO3 (:,K) 
!     ENDDO
      ENDDO

!   UDO3,URO3 ARE IN UNITS OF CM 
!   CFCO2,CFO3 IS THE CONVERSION FACTOR FROM GM/CM2 TO CM

      DO I=1,IMAX
         UR   (I,LP1)=UD   (I,LP1) 
         URCO2(I,LP1)=UDCO2(I,LP1) 
         URO3 (I,LP1)=UDO3 (I,LP1) 
      ENDDO

      DO K=LMAX,1,-1 
      DO I=1,IMAX
         UR   (I,K)=UR   (I,K+1)+VTROW1(I,K)*DIFFCTR
         URCO2(I,K)=URCO2(I,K+1)+VTROW2(I,K)*DIFFCTR
         URO3 (I,K)=URO3 (I,K+1)+VTROW3(I,K)*REFLO3 
      ENDDO
      ENDDO

!   CALCULATE WATER VAPOR TRANSMISSION FUNCTIONS FOR BANDS 2-9;
!   T.F. FOR BAND 1= T.F FOR BAND 2

      DO N=2,NB 
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=ABCFF(N)*UD(I,K) 
         TTU(I,K,N)=ABCFF(N)*UR(I,K) 
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         IF (TTD(I,K,N).GE.50.)  TTD(I,K,N)=50.
         IF (TTU(I,K,N).GE.50.)  TTU(I,K,N)=50.
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=-1.0*TTD(I,K,N)
         TTU(I,K,N)=-1.0*TTU(I,K,N)
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=EXP(TTD(I,K,N))
         TTU(I,K,N)=EXP(TTU(I,K,N))
      ENDDO
      ENDDO
      ENDDO

!   CALCULATE CO2 ABSORPTIONS . THEY WILL BE USED IN BANDS 2-9. 
!   SINCE THESE OCCUPY 50 PERCENT OF THE SOLAR SPECTRUM THE 
!   ABSORPTIONS WILL BE MULTIPLIED BY 2 

      DO I=1,IMAX
         ADCO2(I,1)=0.00
      ENDDO

      DO K=2,LP1
      DO I=1,IMAX
         ADCO2(I,K)=UDCO2(I,K)+0.0129
         ADCO2(I,K)=0.26*LOG(ADCO2(I,K))
         ADCO2(I,K)=EXP(ADCO2(I,K))
         ADCO2(I,K)=2.35E-3*ADCO2(I,K)-7.58265E-4 
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         AUCO2(I,K)=URCO2(I,K)+0.0129
         AUCO2(I,K)=0.26*LOG(AUCO2(I,K))
         AUCO2(I,K)=EXP(AUCO2(I,K))
         AUCO2(I,K)=2.35E-3*AUCO2(I,K)-7.58265E-4 
      ENDDO
      ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        ACO2(I,K)=UCO2(I,K)+0.0129
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        ACO2(I,K)=LOG(ACO2(I,K))
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        ACO2(I,K)=0.26*ACO2(I,K)
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        ACO2(I,K)=EXP(ACO2(I,K))
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        ACO2(I,K)=2.35E-3*ACO2(I,K)-7.58265E-4 
!     ENDDO
!     ENDDO

      DO I=1,IMAX
         AUCO2(I,LP1)=ADCO2(I,LP1) 
      ENDDO

      DO K=1,LP1 
      DO I=1,IMAX 
         ADCO2(I,K)=2.0*ADCO2(I,K) 
         AUCO2(I,K)=2.0*AUCO2(I,K) 
      ENDDO
      ENDDO

!   NOW CALCULATE OZONE ABSORPTIONS. THESE WILL BE USED IN
!   BAND 1. AS THIS OCCUPIES 50 PERCENT OF THE SOLAR SPECTRUM 
!   THE OZONE ABSORPTIONS WILL BE MULTIPLIED BY 2

      DO I=1,IMAX
         ABSDO3(I,1)=0.00 
      ENDDO

      H103P6=103.6*103.6*103.6

      DO K=2,LP1
      DO I=1,IMAX
         ABSDO3(I,K)=1.0+138.6*UDO3(I,K) 
         ABSDO3(I,K)=-0.805*LOG(ABSDO3(I,K))
         ABSDO3(I,K)=EXP(ABSDO3(I,K))
         ABSDO3(I,K)=1.082*UDO3(I,K)*ABSDO3(I,K)+ &
            0.0658*UDO3(I,K)/(1.0+H103P6*UDO3(I,K)*UDO3(I,K)*UDO3(I,K))+ &
            0.02118*UDO3(I,K)/(1.0+0.042*UDO3(I,K)+ &
            0.000323*UDO3(I,K)*UDO3(I,K))
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         ABSUO3(I,K)=1.0+138.6*URO3(I,K) 
         ABSUO3(I,K)=-0.805*LOG(ABSUO3(I,K))
         ABSUO3(I,K)=EXP(ABSUO3(I,K))
         ABSUO3(I,K)=1.082*URO3(I,K)*ABSUO3(I,K)+ &
            0.0658*URO3(I,K)/(1.0+H103P6*URO3(I,K)*URO3(I,K)*URO3(I,K))+ &
            0.02118*URO3(I,K)/(1.0+0.042*URO3(I,K)+ &
            0.000323*URO3(I,K)*URO3(I,K))
      ENDDO
      ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        AO3(I,K)=1.0+138.6*UO3(I,K) 
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        AO3(I,K)=LOG(AO3(I,K))
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        AO3(I,K)=-0.805*AO3(I,K)
!     ENDDO
!     ENDDO

!     DO K=2,LLP1
!     DO I=1,IMAX
!        AO3(I,K)=EXP(AO3(I,K)) 
!     ENDDO
!     ENDDO

!     H103P6=103.6*103.6*103.6
!     DO K=2,LLP1
!     DO I=1,IMAX 
!        AO3(I,K)=1.082*UO3(I,K)*AO3(I,K)+ 
!    &        0.0658*UO3(I,K)/(1.0+H103P6*UO3(I,K)*UO3(I,K)*UO3(I,K))+
!    &        0.02118*UO3(I,K)/(1.0+0.042*UO3(I,K)+ 
!    &        0.000323*UO3(I,K)*UO3(I,K))
!     ENDDO
!     ENDDO

      DO I=1,IMAX
         ABSUO3(I,LP1)=ABSDO3(I,LP1) 
      ENDDO

!     DO K=1,LLP2
      DO K=1,LP1
      DO I=1,IMAX 
!        AO3(I,K)=2.0*AO3(I,K)
         ABSDO3(I,K)=2.0*ABSDO3(I,K) 
         ABSUO3(I,K)=2.0*ABSUO3(I,K) 
      ENDDO
      ENDDO

!     WRITE (*,101) ((K,UD(IP,K),UR(IP,K),K=1,LP1),IP=1,IMAX) 
!     WRITE (*,105) ((K,UDO3(IP,K),URO3(IP,K),ABSDO3(IP,K), 
!    1 ABSUO3(IP,K),K=1,LP1),IP=1,IMAX) 
!     WRITE (*,105) ((K,UDCO2(IP,K),URCO2(IP,K),ADCO2(IP,K),
!    1 AUCO2(IP,K),K=1,LP1),IP=1,IMAX)

!   COMBINE ABSORPTIONS AND TRANSMISSIONS TO OBTAIN A 
!   TRANSMISSION FUNCTION FOR EACH OF THE 9 BANDS.

      DO K=1,LP1 
      DO I=1,IMAX 
         TTD(I,K,1)=TTD(I,K,2)*(1.0-ABSDO3(I,K))
      ENDDO
      ENDDO

      DO K=1,LMAX 
      DO I=1,IMAX 
         TTU(I,K,1)=TTU(I,K,2)*(1.0-ABSUO3(I,K))
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1 
      DO I=1,IMAX 
         TTD(I,K,N)=TTD(I,K,N)*(1.0-ADCO2(I,K)) 
      ENDDO
      ENDDO
      ENDDO

      DO N=1,NB
      DO I=1,IMAX 
         TTU(I,LP1,N)=TTD(I,LP1,N) 
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LMAX 
      DO I=1,IMAX 
         TTU(I,K,N)=TTU(I,K,N)*(1.0-AUCO2(I,K)) 
      ENDDO
      ENDDO
      ENDDO

!   IN THE 850 LOOP BELOW AND THE 855 LOOP,WE WILL SCALE DFN(IP,1,N)
!   TO UNITY. AFTER THE 850 LOOP,WE MULTIPLY BY DFN(IP,1,N) FROM THE
!   6 LOOP TO GET THE ACTUAL DFN'S AND UFN'S.
!   THE 855 LOOP=SCALED DOWNWARD FLUX ADOVE TOPMOST CLOUD(REDUN-
!   DANTLY OBTAINED TO THE GROUND) .ALSO,UFN IS INITIALIZED TO 0

      DO N=1,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFN(I,K,N)=TTD(I,K,N)
         UFN(I,K,N)=0.00
      ENDDO
      ENDDO
      ENDDO
 
!***EVALUATE LOGICAL ARRAYS USED FOR CLOUD CALCULATIONS
!----BCLDS : TRUE IF KK IS < OR = TO NCLDS (NO. CLOUDS AT GRID PT I)
!----BJTOP : TRUE IF K IS < OR = TO KTOPSW(I,2) (INDEX OF TOP CLOUD,
!                    IF ANY; LP1 IS NO CLOUD AT GRID PT)
 
      DO KK=1,LP1
      DO I=1,IMAX
         BCLDS(I,KK)=KK.LE.NCLDS(I)
      ENDDO
      ENDDO

      DO K=1,LP1
      DO I=1,IMAX
         BJTOP(I,K)=K.LE.KTOPSW(I,2)
      ENDDO
      ENDDO

!---COUNT NO. OF PTS IN EACH ROW FOR WHICH BCLDS IS TRUE (ICNT)
!   AND FOR WHICH BJTOP IS TRUE (ICNT1)
 
      DO K=1,LP1
         ICNT(K)=0
         ICNT1(K)=0
         DO I=1,IMAX
            IF (BCLDS(I,K)) THEN
               ICNT(K)=ICNT(K)+1
            ENDIF
            IF (BJTOP(I,K)) THEN
               ICNT1(K)=ICNT1(K)+1
            ENDIF
         ENDDO
      ENDDO

!---FIND NO. OF CLOUD LEVELS WITH NONZERO VALUES OF ICNT
!---FIND NO. OF PRESSURE LEVELS WITH NONZERO VALUES OF ICNT1
 
      KCLDS=0
      KJTOP=0
      DO K=1,LP1
         IF (ICNT(K).GT.0) THEN
            KCLDS=KCLDS+1 
         ENDIF
         IF (ICNT1(K).GT.0) THEN
            KJTOP=KJTOP+1
         ENDIF
      ENDDO

!***IF NO CLOUDS AT ALL EXIST IN THE ROW, THE CALCULATIONS ARE
!   DRASTICALLY SIMPLIFIED.

!MPP  IF (KCLDS.EQ.0) THEN
      IF (KCLDS.EQ.-1) THEN
         DO N=1,NB

            IF (N.EQ.1) THEN
               DO I=1,IMAX
                  REFL(I)=CUVRF(I,2)
!CCCC             REFL(I)=SALB(I)
                  RRAY(I)=0.219/(1.0+0.816*COSZRO(I))
                  REFL(I)=RRAY(I)+(1.0-RRAY(I))*(1.0-RRAYAV)*REFL(I)/ &
                          (1.0-REFL(I)*RRAYAV)
                  ALFA(I)=REFL(I)
               ENDDO
            ELSE
               DO I=1,IMAX
                  ALFA(I)=CIRRF(I,2)
!CCCC             ALFA(I)=SALB(I)
               ENDDO
            ENDIF

            DO I=1,IMAX
               VV(I)=ALFA(I)*DFN(I,LP1,N)/TTU(I,LP1,N)
            ENDDO

            DO K=1,LP1
            DO I=1,IMAX
               UFN(I,K,N)=VV(I)*TTU(I,K,N)
            ENDDO
            ENDDO

         ENDDO
      ENDIF

!***********************************************************************
!     ****** COMPUTE NORMAL CASE: AT LEAST 1 PT HAS A CLOUD ******

!MPP                  IF (KCLDS.NE.0) THEN
                      IF (KCLDS.GE.0) THEN

!***********************************************************************

!---FIND  HIGHEST PRESSURE LEVEL WITH AT LEAST 1 PT BELOW TOP CLOUD
      KCLDS2=0
      DO K=1,LP1
         IF (ICNT1(K).EQ.IMAX) THEN
            KCLDS2=KCLDS2+1
         ENDIF
      ENDDO
      KCLDS2=KCLDS2+1

!    -------------------
      DO 2105 KK=1,KCLDS
!    -------------------
!---DETERMINE WHETHER A CLOUD LAYER KK HAS AT LEAST 1 GRID PT WITH A
!   "THICK CLOUD"
!---DETERMINE WHETHER THERE IS AT LEAST 1 GRID PT WHERE THE BOTTOM
!  OF CLOUD LAYER KK DOES NOT COINCIDE WITH THE TOP OF CLOUD LAYER
!  KK+1

      IINCL(KK)=0
      IBETCL(KK)=0
      DO K=KCLDS2,LP1
      DO I=1,IMAX
         IF (BCLDS(I,KK) .AND. &
           (K.GT.KTOPSW(I,KK+1) .AND. K.LT.KBTMSW(I,KK+1))) THEN
               IINCL(KK)=IINCL(KK)+1
         ENDIF
         IF (BCLDS(I,KK) .AND. &
           (K.GE.KBTMSW(I,KK+1) .AND. K.LE.KTOPSW(I,KK+2))) THEN
               IBETCL(KK)=IBETCL(KK)+1
         ENDIF
      ENDDO
      ENDDO
!    -------------------
2105  CONTINUE
!    -------------------

!***COMPUTE VISIBLE BAND GROUND REFLECTIVITY USING LACIS-HANSEN 
!   PARAMETERIZATION

      DO I=1,IMAX
         REFL(I)=CUVRF(I,NCLDS(I)+2)
!CCCC    REFL(I)=SALB(I)
      ENDDO

      DO IP=1,IMAX
         RRAY(IP)=0.219/(1.0+0.816*COSZRO(IP))
         REFL(IP)=RRAY(IP)+(1.0-RRAY(IP))*(1.0-RRAYAV)*REFL(IP)/ &
                           (1.0-REFL(IP)*RRAYAV)
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         PPTOP(I,KK)=PP(I,KTOPSW(I,KK+1)) 
         PPBOT(I,KK)=PP(I,KBTMSW(I,KK+1))
      ENDDO
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         IF (PPTOP(I,KK).NE.PPBOT(I,KK)) THEN
            DPCLD(I,KK)=1.0/(PPTOP(I,KK)-PPBOT(I,KK))
         ELSE
            DPCLD(I,KK)=0.00
         ENDIF
      ENDDO
      ENDDO

!***WE NOW OBTAIN AN INDEX FOR (I,NCLDS(I)+1-KK).WE FORCE THIS
!   INDEX TO HAVE A MINIMUM VALUE IN THE (0,IMAX) RANGE.

      DO KK=1,KCLDS+1
      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            INDXK(I)=KK
         ELSE
            INDXK(I)=NCLDS(I)
         ENDIF
         INDX4(I,KK)=(NCLDS(I)-INDXK(I))*IMAX+I
      ENDDO
      ENDDO

!-----------------------------------------------------------------------
!   THE REST OF THE CLOUD CALCULATION IS PERFORMED INSIDE A
!   BAND (FREQUENCY) LOOP OVER N, RUNNING FROM 1 TO NB
!-----------------------------------------------------------------------

                         DO 2301 N=1,NB

!     print *, 'NB,N=',NB,N
!-----------------------------------------------------------------------

!***INITIALIZE CR TO ZERO AND CT TO ONE***
      DO K=1,LP1
      DO I=1,IMAX
         CR(I,K)=0.00
         CT(I,K)=1.0
      ENDDO
      ENDDO

!***OBTAIN CLOUD REFLECTION AND TRANSMISSION COEFFICIENTS, FIRST FOR
!   VISIBLE BAND (N=1) THEN FOR NEAR IR BANDS (N=2-NB) 
!---FIRST, THE VISIBLE BAND:

!               ----------------
                IF (N == 1) THEN
!               ----------------

      DO KK=1,KCLDS
      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            CR(I,KK+1)=CUVRF(I,KK+1)*CAMT(I,KK+1)
            CT(I,KK+1)=1.0-CR(I,KK+1)
         ENDIF
      ENDDO
      ENDDO

!---USE THIS INDEX FOR SPECIFYING VISIBLE BAND GROUND ALBEDO 
!   AS REFL(I):

      DO I=1,IMAX
         CR(I,NCLDS(I)+2)=REFL(I)
      ENDDO

!               ----------------
                     ENDIF
!               ----------------

!---NOW, THE NEAR IR BANDS; HERE THE GROUND CAN BE HANDLED AS PART
!   OF THE CLOUD LOOP

!               ----------------
                IF (N > 1) THEN
!               ----------------

      DO I=1,IMAX 
          CR(I,2)=CIRRF(I,2)*CAMT(I,2)
          CT(I,2)=1.0-CAMT(I,2)*(CIRRF(I,2)+CIRAB(I,2))
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
!       print *, 'I,KK,BCLDS,CIRRF,CAMT,CIRAB=',I,KK,
!    &           BCLDS(I,KK),CIRRF(I,KK+2),CAMT(I,KK+2),CIRAB(I,KK+2)
         IF (BCLDS(I,KK)) THEN
            CR(I,KK+2)=CIRRF(I,KK+2)*CAMT(I,KK+2)
            CT(I,KK+2)=1.0-CAMT(I,KK+2)*(CIRRF(I,KK+2)+CIRAB(I,KK+2))
         ENDIF
      ENDDO
      ENDDO

!               ----------------
                     ENDIF
!               ----------------

      DO K=1,LP1
      DO I=1,IMAX
         ALFAU(I,K)=0.00
      ENDDO
      ENDDO


!***FOR EXECUTION OF THE CLOUD LOOP, IT IS NECESSARY TO SEPARATE OUT
!   TRANSMISSION FCTNS AT THE TOP AND BOTTOM OF THE CLOUDS, FOR
!   EACH BAND N. THE REQUIRED QUANTITIES ARE:
!      TTD(I,KTOPSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+2: 
!      TTD(I,KBTMSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+1: 
!      TTU(I,KTOPSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+2:
!      AND INVERSES OF THE ABOVE. THE ABOVE QUANTITIES ARE STORED 
!      IN TDCL1,TDCL2,TUCL1,TDCL1I,TDCL2I,TUCLI,RESPECTIVELY, AS
!      THEY HAVE MULTIPLE USE IN THE PGM.

!---COMPUTE GATHERS
      DO KK=1,KCLDS+1
      DO I=1,IMAX
         TDCL1(I,KK)=TTD(I,KTOPSW(I,KK+1),N)
         TUCL1(I,KK)=TTU(I,KTOPSW(I,KK+1),N)
      ENDDO
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         TDCL2(I,KK)=TTD(I,KBTMSW(I,KK+1),N)
      ENDDO
      ENDDO

!---COMPUTE INVERSES
      DO KK=1,KCLDS
      DO I=1,IMAX
         TDCL2I(I,KK)=1.0/TDCL2(I,KK)
      ENDDO
      ENDDO

      DO KK=1,KCLDS+1
      DO I=1,IMAX
        TDCL1I(I,KK)=1.0/TDCL1(I,KK)
        TUCL1I(I,KK)=1.0/TUCL1(I,KK)
      ENDDO
      ENDDO


!   TCLU(LL) IS TRANSMISSION FUNCTION FROM TOP OF NEXT LOWER
!   CLOUD TO TOP OF UPPER CLOUD. TCLD(LL) IS T.F. FROM TOP
!   OF NEXT LOWER CLOUD TO BOTTOM OF UPPER CLOUD. LL=NC1 IS
!   THE LOWEST BOTTOM CLOUD (THE GROUND) ; LL=1 IS THE
!   HIGHEST UPPER CLOUD. 

         TCLU = 0.0
         TCLD = 0.0
      DO KK=1,KCLDS
      DO I=1,IMAX
         TCLU(I,KK+1)=TDCL1(I,KK+1)*TDCL1I(I,KK)*CT(I,KK+1)
         TCLD(I,KK+1)=TDCL1(I,KK+1)*TDCL2I(I,KK) 
      ENDDO
      ENDDO

!***WE DEFINE TCLD (I,1) AS TTD(I,KTOPSW(I,2),N)
      DO I=1,IMAX
         TCLD(I,1)=TDCL1(I,1)
      ENDDO

      DO I=1,IMAX
         DFNCLU(I,1)=TCLD(I,1)
      ENDDO


!   THE FOLLOWING CALCULATION IS FOR ALFAT: THE RATIO BETWEEN
!   THE DOWNWARD FLUX AT THE TOP OF THE HIGHEST CLOUD (IF
!   PRESENT) OR THE GROUND TO THE UPWARD FLUX AT THE SAME
!   LEVEL, TAKING INTO ACCOUNT MULTIPLE REFLECTIONS FROM
!   CLOUDS, IF PRESENT

!  --- Reshape 2-D arrays to 1-D ---
      DO K=1,LP1
      DO I=1,IMAX
        I1=(K-1)*IMAX+I
        CR1D(I1)=CR(I,K)
        ALFU1D(I1)=ALFAU(I,K)
        TCLD1D(I1)=TCLD(I,K)
        TCLU1D(I1)=TCLU(I,K)
      ENDDO
      ENDDO

!     print *, 'KCLDS=',KCLDS

      DO KK=1,KCLDS
!     -------------

      DO I=1,IMAX
         TEMP1(I)=CR1D(INDX4(I,KK)+2*IMAX)
         TEMP2(I)=ALFU1D(INDX4(I,KK)+IMAX)
         TEMP3(I)=TCLU1D(INDX4(I,KK)+IMAX)
         TEMP4(I)=CR1D(INDX4(I,KK)+IMAX)
         TEMP5(I)=TCLD1D(INDX4(I,KK)+IMAX)
      ENDDO
!     print *, 'TEMP1=',TEMP1

      DO I=1,IMAX
         TEMP6(I)=(TEMP1(I)+TEMP2(I))*TEMP3(I)*TEMP3(I) / &
                  (1.0-(TEMP1(I)+TEMP2(I))*TEMP4(I)*TEMP5(I)*TEMP5(I))
      ENDDO
!     print *, 'TEMP6=',TEMP6

      DO I=1,IMAX
         ALFU1D(INDX4(I,KK))=TEMP6(I)
      ENDDO

!  --- Reshape 1-D array into 2-D array -----
      DO K=1,LP1
      DO I=1,IMAX
         I1=(K-1)*IMAX+I
         ALFAU(I,K)=ALFU1D(I1)
      ENDDO
      ENDDO

!     print *, 'ALFAU=',ALFAU
!     -------------
      ENDDO

!***DEFINE ALFA FROM ALFAU(I,1) AND CR(I,2):
!***ALFA IS THE SYSTEM REFLECTION COEFFICIENT ABOVE THE TOP CLOUD
!    (OR GROUND, IF NO CLOUD AT GRID PT I )

      DO I=1,IMAX
         ALFA(I)=ALFAU(I,1)+CR(I,2)
      ENDDO

!     UPWARD FLUX ABOVE TOPMOST CLOUD
      DO I=1,IMAX
         UFNCLU(I,1)=TCLD(I,1)*ALFA(I)
      ENDDO

      DO I=1,IMAX
         TEMP2(I)=TUCL1I(I,1)*UFNCLU(I,1)
      ENDDO
!     print *, 'TEMP2=',TEMP2

      DO K=1,KJTOP      
      DO I=1,IMAX
         IF (BJTOP(I,K)) THEN
            UFN(I,K,N)=TEMP2(I)*TTU(I,K,N)
         ENDIF
      ENDDO
      ENDDO
!     print *, 'UFN=',UFN

!   CALCULATE UFN AT CLOUD TOPS AND DFN AT CLOUD BOTTOMS

      DO I=1,IMAX
         VV(I)=1.0
      ENDDO

      DO KK=1,KCLDS 
!     -------------
!     print *, 'KK,KCLDS=',KK,KCLDS
!     print *, 'TCLU=',TCLU

      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            UFNCLU(I,KK+1)=ALFAU(I,KK)*VV(I)*TCLD(I,KK)/TCLU(I,KK+1)
            DFNCLD(I,KK)=VV(I)*TCLD(I,KK)* &
                         TCLU(I,KK+1)*TDCL2(I,KK)*TDCL1I(I,KK+1) + &
                         UFNCLU(I,KK+1)*TCLD(I,KK+1)*CR(I,KK+1)
         ELSE
            UFNCLU(I,KK+1)=UFN(I,LP1,N)
            DFNCLD(I,KK)=DFN(I,LP1,N)
         ENDIF
      ENDDO
!     print *, 'UFNCLU=',UFNCLU
!     print *, 'DFNCLD=',DFNCLD

      DO I=1,IMAX
         VV(I)=DFNCLD(I,KK)
      ENDDO

!     print *, 'VV=',VV
!     print *, 'KTOPSW(KK+2)=',KTOPSW(:,KK+2)
!     print *, 'KBTMSW(KK+1)=',KBTMSW(:,KK+1)

      DO I=1,IMAX
         UFN(I,KTOPSW(I,KK+2),N)=UFNCLU(I,KK+1)
         DFN(I,KBTMSW(I,KK+1),N)=DFNCLD(I,KK)
      ENDDO

!     -------------
      ENDDO


!     NOW OBTAIN DFN AND UFN FOR LEVELS BETWEEN THE CLOUDS
      DO 2401 KK=1,KCLDS
!---SKIP IF THERE ARE NO SPACES BETWEEN CLOUD LAYERS KK AND KK+1,
!    FOR ANY GRID PT:
      IF (IBETCL(KK).EQ.0) GO TO 2401
         DO K=KCLDS2,LP1
         DO I=1,IMAX
            IF (BCLDS(I,KK) .AND. &
               (K.GE.KBTMSW(I,KK+1) .AND. K.LE.KTOPSW(I,KK+2))) THEN
                  UFN(I,K,N)=UFNCLU(I,KK+1)*TTU(I,K,N)*TUCL1I(I,KK+1)
                  DFN(I,K,N)=DFNCLD(I,KK)*TTD(I,K,N)*TDCL2I(I,KK)
            ENDIF
         ENDDO
         ENDDO
2401  CONTINUE


!     NOW OBTAIN DOWNWARD AND UPWARD FLUXES FOR LEVELS,IF ANY,
!     BETWEEN THE TOPS AND BOTTOMS OF CLOUDS. THE ASSUMPTION OF
!     CONSTANT HEATING RATE IN THESE REGIONS IS USED.

!***OBTAIN FLUXES AT TOP AND BOTTOM OF CLOUDS
      DO 2501 KK=1,KCLDS
!---SKIP IF THERE ARE NO "THICK CLOUDS" AT ALL IN CLOUD LEVEL KK
      IF (IINCL(KK).EQ.0) GO TO 2501
!
!***OBTAIN DOWNWARD FLUXES AT CLOUD TOPS AND UPWARD FLUXES AT
!   CLOUD BOTTOMS

      IF (KK.GT.1) THEN
         DO I=1,IMAX
            DFNCLU(I,KK)=DFN(I,KTOPSW(I,KK+1),N)
         ENDDO
      ENDIF

      DO I=1,IMAX
         UFNCLD(I,KK)=UFN(I,KBTMSW(I,KK+1),N)
      ENDDO

      DO I=1,IMAX
         TEMP1(I)=(UFNCLU(I,KK)-UFNCLD(I,KK))*DPCLD(I,KK)
         TEMP2(I)=(DFNCLU(I,KK)-DFNCLD(I,KK))*DPCLD(I,KK)
      ENDDO


      DO K=KCLDS2,LP1
      DO I=1,IMAX
         IF (BCLDS(I,KK) .AND. &
             (K.GT.KTOPSW(I,KK+1) .AND. K.LT.KBTMSW(I,KK+1))) THEN
                UFN(I,K,N)=UFNCLU(I,KK)+TEMP1(I)*(PP(I,K)-PPTOP(I,KK))
                DFN(I,K,N)=DFNCLU(I,KK)+TEMP2(I)*(PP(I,K)-PPTOP(I,KK))
         ENDIF
      ENDDO
      ENDDO

2501  CONTINUE

!-----------------------------------------------------------------------

2301                           CONTINUE

!-----------------------------------------------------------------------

                                ENDIF

!***********************************************************************

!   CALCULATE ENTERING FLUX AT THE TOP
!   LOOP 860 SCALES THE DFN'S AND UFN'S TO THE CORRECT DFN(I,1,N)

      DO N=1,NB
      DO I=1,IMAX
         DFNTOP(I,N)=SSOLAR(I)*6.97667E5*PWTS(N)
!OLD     DFNTOP(I,N)=SSOLAR*6.97667E5*COSZRO(I)*TAUDAR(I)*PWTS(N)
      ENDDO
      ENDDO

      DO N=1,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFN(I,K,N)=DFN(I,K,N)*DFNTOP(I,N)
         UFN(I,K,N)=UFN(I,K,N)*DFNTOP(I,N)
      ENDDO
      ENDDO
      ENDDO

!   SUM OVER BANDS

      DO K=1,LP1
      DO I=1,IMAX
         DFSW(I,K)=DFN(I,K,1)
         UFSW(I,K)=UFN(I,K,1)
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFSW(I,K)=DFSW(I,K)+DFN(I,K,N)
         UFSW(I,K)=UFSW(I,K)+UFN(I,K,N)
      ENDDO
      ENDDO
      ENDDO

      DO K=1,LP1
      DO I=1,IMAX
         FSW(I,K)=UFSW(I,K)-DFSW(I,K)
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         HSW(I,K)=RADCON*(FSW(I,K+1)-FSW(I,K))/DP(I,K)
      ENDDO
      ENDDO

!-----------------------------------------------------------------------

      END SUBROUTINE SWRAD_ORIG

!#######################################################################
      SUBROUTINE SHORTWAVE_INIT

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      END SUBROUTINE SHORTWAVE_INIT

!#######################################################################

      SUBROUTINE SHORTWAVE_END

      module_is_initialized = .false.

!---------------------------------------------------------------------

      END SUBROUTINE SHORTWAVE_END

!#######################################################################

      END MODULE SHORTWAVE_MOD
