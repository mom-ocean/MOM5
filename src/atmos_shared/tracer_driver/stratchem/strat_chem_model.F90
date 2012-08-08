MODULE STRAT_CHEM_MOD
! automatic conversion to free f90 compatible form 
! free.pl strat_chem_model.f  
! linewidth: 72
! file names: strat_chem_model.f 
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
!/
!/ THIS MODIFICATION SET FOR THE UNIFIED MODEL IS A 
!/ REAL STRATOSPHERIC CHEMICAL MODEL         
!/ WRITTEN BY J AUSTIN (10/02/94); MODIFIED 26/7/95 FOR 49 LEVEL MODEL  
!/ REVISED FOR CLIMATE LENGTH INTEGRATIONS 24/9/96
!/ REVISED FOR PARALLEL CODE 13/11/97
!/ REVISED CHEMICAL MODEL 14/1/98
!/ REVISED CHEMICAL MODEL WITH MAINZ ANALYTIC SCHEME 1/5/99
!/ REVISED INTEGRATION SCHEME AFTER STEILET AL. 1998 14/6/99
!/ UPDATED TO VN4.5 64 LEVELS J. AUSTIN 27/1/00
!/ UPDATED CH4 OXIDATION, TROPOSPHERIC CHEMISTRY 14/2/02
!/ UPDATED FOR USE IN THE FMS, 15/12/03
!/ REVISED, 13/10/04: SIMPLIFIED PSCs, IMPROVED LONG-LIVED TRACERS,
!/      AGE OF AIR PARAMETERISATION FOR HALOGEN RATES OF CHANGE
!/

implicit none
private

public CHEMISTRY, ZEN2, DCLY_DT, SEDIMENT

contains

      SUBROUTINE CHEMISTRY (ALON, ALAT, JLZ, KL, DY, H2O, DAGESQ, OZCOL, &
        PRESS, RHO, TEMP, VTEMP, COZEN, CDX, CHLB, OZB, BNAT, BICE, PHOTO, &
        SOLARDATA, OZON, COSP, COSPHC, ANOY, AEROSOL, DT, MERIDS, CH_TEND, &
        OZONE, O3_PROD, H2O_TEND, MYPE, ITIME)
! """" ----------------------------------------------------------------*
! VERSION NUMBER   : 1   DATE:  6/10/1989.
!   NOTES          : CODE REWRITTEN TO SOLVE EQUATIONS MORE DIRECTLY
!                  : IN FAMILY MODE TO COMPARE MODEL RESULTS
! VERSION NUMBER   : 2   DATE:  1/06/1990.
!   NOTES          : REACTION SET SIMPLIFIED TO REDUCE CPU TIME
! VERSION NUMBER   : 3   DATE: 02/04/1993.
!   NOTES          : REVISED HETEROGENEOUS REACTION SCHEME
! VERSION NUMBER   : 4   DATE: 09/02/1994ET SEQ.
!   NOTES          : ADAPTED FOR USE IN THE UNIFIED MODEL
! VERSION NUMBER   : 5   DATE: 15/12/2003
!   NOTES          : ADAPTED FOR USE IN THE FMS               
! VERSION NUMBER   : 6   DATE: 23/03/2006
!   NOTES          : CONVERTED TO f90 FOR USE IN THE FMS               
! """" ----------------------------------------------------------------*
! PURPOSE          : TO INTEGRATE A ZERO-DIMENSION(POINT) CHEMICAL      
!                    MODEL ALONG A PREDETERMINED TRAJECTORY FOR A       
!                    TIME PERIOD IDT
! METHOD           : USES THE EULER BACKWARD INTEGRATION SCHEME         
! INPUT/OUTPUT LIST:
!     JLX          = LAT.  NUMBER (NEEDED FORERROR DIAGNOSTICS)        
!     KL           = LEVEL NUMBER (NEEDED FOR PHOTOLYSIS RATES)         
!     DT           = INTEGRATION TIMESTEP
!     DY           = ARRAY OF FAMILY CONCENTRATIONS
!     OZCOL        = OZONE COLUMN ABOVE PARTICLE
!     FACT         = RATIO OF OZONE TO OY (FOR PHOTOLYSIS RATES)        
!     RHO          = AIR DENSITY AT THE PARTICLE
!     TEMP         = TEMPERATURE AT THE PARTICLE
!     CHEMINC      = INCREMENTS TO FAMILY CONCENTRATIONS      
!                    AFTER THIS CALL TO OAERO
!
!  SUBPROGRAM CALLS:
!     REACTION_RATES  - CALCULATES REACTION RATES (CALLED ONCE) 
!     PSC/GAMMA_AER/HETRATES - CALCULATE HETEROGENEOUS REACTION RATES
!     PHOTO   - CALCULATES PHOTODISSOCIATION RATES (CALLED ONCE)        
!----------------------------------------------------------------------*
      IMPLICIT NONE
      

      integer ::  nchem, jimpl, ihet, nhet, irx, irx4
      PARAMETER (NCHEM=21,JIMPL=18,IHET=72,NHET=9)
      PARAMETER (IRX=33,IRX4=IRX*4) 

      integer, intent(in) :: MERIDS, JLZ, KL, MYPE
      integer, intent(in) :: ITIME(6)
      real, intent(in)    :: CDX, DT
      real, intent(in), dimension(MERIDS) :: ALON, ALAT,   &
                           DAGESQ, OZCOL, PRESS, RHO, &
                           TEMP, VTEMP, COZEN, &
                           AEROSOL
      real, intent(in), dimension(MERIDS,NCHEM) :: DY
      real, intent(in), dimension(90,15) :: CHLB
      real, intent(in), dimension(144,90) :: OZB
      real, intent(in) :: PHOTO(132,14,11,48)
      real, intent(in) :: SOLARDATA(1801)
      real, intent(in) :: OZON(11,48)
      real, intent(inout) :: COSP(14)
      real, intent(inout) :: COSPHC(48)
      real, intent(in) :: ANOY(90,48)

      real, intent(inout), dimension(MERIDS)       :: H2O, BNAT, BICE, OZONE, O3_PROD, H2O_TEND
      real, intent(inout), dimension(MERIDS,NCHEM) :: CH_TEND
              
!----------------------------------------------------------------------*
! VARIABLES PASSED TO OR FROM HIGHER SUBROUTINES
!----------------------------------------------------------------------*


      INTEGER :: ISTEP_CHEM, IC, IL, ILX, JLX, IZ, NC, KLX, IR, &
                 ITX, ITNO, ITIME_LEFT, IT, IM, IDX, ICYC, IDAYS, LEAPD,&
                 JMOD, ISC
      REAL :: AGESQ
      INTEGER :: IMON(12)
!----------------------------------------------------------------------*
! VARIABLES FOR THIS SUBROUTINE
!----------------------------------------------------------------------*
#define _ALLOC

#ifdef _ALLOC
     REAL, ALLOCATABLE, DIMENSION(:,:) :: RATES, FOTO
     REAL, ALLOCATABLE, DIMENSION(:,:) :: DYY, DYZ, DQT, DPT, RAT, &
              DYZ1, DYY0, DYY1, COND, GAMMA
#else
     REAL :: RATES(MERIDS,IHET+NHET-1)
     REAL :: FOTO(MERIDS,IRX4)
     REAL :: DYY(MERIDS,NCHEM)
     REAL :: DYZ(MERIDS,JIMPL)
     REAL :: DQT(MERIDS,NCHEM)
     REAL :: DPT(MERIDS,NCHEM)
     REAL :: RAT(MERIDS,3)
     REAL :: DYZ1(MERIDS,JIMPL)
     REAL :: DYY0(MERIDS,NCHEM)
     REAL :: DYY1(MERIDS,NCHEM)
     REAL :: COND(MERIDS,3)
     REAL :: GAMMA(MERIDS,NHET)
#endif
     REAL, DIMENSION(MERIDS) ::  CLOX, BROX, ANOZ, H2O1, H2O0, DP1, DQ1, &
            H2SO4, TICE, WH2SO4, AM, AW, ALIQ, RMEAN, ASAT, RNAT, RICE, RTT

     REAL :: DVDT, CH4, AN2O, CH4TR, AN2OTR, ANAT, AICE, VTAU10, VTAU100, &
            VTAUYR, FACT, DELT, SUM, RATIO1, RATIO2, ANUM, ADEN, ANOX, DAX, &
            DA1, DA2, DD1, DB1, DB2, DD2, DDNUM, DA3, DC2, DC3, DD3, R1, R2, &
            DD13, VDDNUM, DAC2, DA12, DD12, CLX, TAU, RAT1, BX1, AX, BX, CX, &
            XYZ, DQT1, DQT2, DQT3, T0, T1, BRX, DIFF, VRHO, BX2, &
            RATIO, X1, X2, SOL, FRAC, SOLARTIME, SOLAR27, DX1, D11, D27, F11, &
            F27, PI, O2, H2, minusT0
      DATA  PI/3.141592653589793/
      DATA  O2/0.2095/,H2/0.5E-6/                                
      DATA  ISTEP_CHEM/900/ ! Chemical timestep (s)
      DATA  ANAT,AICE/10.0, 10.0/
      DATA  IMON/0,31,59,90,120,151,181,212,243,273,304,334/
!
!  ISC IS THE SRES SCENARIO FOR CH4 AND N2O
!  1 = A1B, 2 = A2, 3 = B1, 4 = timeslice
!
      DATA  ISC/1/
! ---------------------------------------------------------------------*
!                                                                       
!   THE ARRAYS DYY AND DYZ CONTAIN THE SPECIES AMOUNTS.  THE SPECIES    
!   DYY ARE CALCULATEDEXPLICITLY, WHILE THE SPECIES DYZ ARE FOUND      
!   IMPLICITLY USING PHOTOCHEMICALEQUILIBRIUM AND RELATED APPROACHES.  
!   THE VALUES IN DYZ ARE RECALCULATED AFRESH FOREVERY CALL TO OAERO.  
!   THE ORDER OF THE SPECIES CONTAINED IN THE ARRAYS DYY AND DYZ IS:    
!      DYY                                                              
!    1. HNO3     2. N2O5   3. H2O2    4. HCL    5. HOCL  6. CLONO2      
!    7. H2CO     8. OY     9. HOBR   10. HNO4  11. HBR  12. BRONO2      
!   13. CH3OOH  14. CO    15. NAT+NY 16. CLY   17. BRY  18. CH4  
!   19. Strat H2O 20. N2O 21. AGE  
!                                                                       
!      DYZ                                                              
!    1.   O3      2.  O3P    3.  O1D    4.  NO      5.  NO2    6. NO3   
!    7.   OH      8.  HO2    9.  H     10.  CL     11.  CLO   12. BR    
!   13.   BRO    14.  CL2O2 15.  BRCL  16.  CH3O2  17.  HONO  18. N
!----------------------------------------------------------------------*

#ifdef _ALLOC
     allocate(RATES(MERIDS,IHET+NHET-1))
     allocate(FOTO(MERIDS,IRX4))
     allocate(DYY(MERIDS,NCHEM))
     allocate(DYZ(MERIDS,JIMPL))
     allocate(DQT(MERIDS,NCHEM))
     allocate(DPT(MERIDS,NCHEM))
     allocate(RAT(MERIDS,3))
     allocate(DYZ1(MERIDS,JIMPL))
     allocate(DYY0(MERIDS,NCHEM))
     allocate(DYY1(MERIDS,NCHEM))
     allocate(COND(MERIDS,3))
     allocate(GAMMA(MERIDS,NHET))
#endif
!
! Set PSCs initially to zero
!
     DO IL = 1,MERIDS
       BNAT(IL) = 0.0
       BICE(IL) = 0.0
     ENDDO
!
! FIND CONCENTRATIONS AND TRENDS OF CH4 AND N2O, ACCORDING TO THE SRES SCENARIO
!
      CALL GHGS(ITIME,ISC,CH4,AN2O,CH4TR,AN2OTR)
!   
!  RELAX TOWARDS SPECIFIC LOWER BOUNDARY VALUES WITH A TIMESCALE OF 10 DAYS
!  (SPECIES 1-15) SPECIES 16-21 ARE TREATED DIFFERENTLY
!                   
      VTAU10 = 1.0/(10.0*86400.0)
      VTAU100 = 0.1*VTAU10
      VTAUYR = 1.0/(365.25*86400.0)
      IF(KL.EQ.48) THEN                                                  
        DO IC = 1,15                                            
          IF(IC.EQ.8) THEN
            DO IL = 1,MERIDS                                            
              ILX = 1 + ALON(IL)*72/PI
              JLX = 1 + (PI*0.5 + ALAT(IL))*89.0/PI
              IF(ILX.GT.144) ILX = ILX - 144 
              CH_TEND(IL,IC) = (OZB(ILX,JLX) - DY(IL,IC))*VTAU10      
            ENDDO
          ELSE             
            DO IL = 1,MERIDS                                            
              JLX = 1 + (PI*0.5 + ALAT(IL))*89.0/PI
              CH_TEND(IL,IC) = (CHLB(JLX,IC) - DY(IL,IC))*VTAU10     
            ENDDO
          ENDIF
        ENDDO                                             
        DO IL = 1,MERIDS
          CH_TEND(IL,16) = - DY(IL,16)*VTAU10     
          CH_TEND(IL,17) = - DY(IL,17)*VTAU10     
          CH_TEND(IL,18) = (CH4 - DY(IL,18))*VTAU10     
          CH_TEND(IL,19) = (3.0E-6 - DY(IL,19))*VTAU10     
          CH_TEND(IL,20) = (AN2O - DY(IL,20))*VTAU10     
        ENDDO
        RETURN
      ENDIF             
!----------------------------------------------------------------------*
!   NO CHANGE TO CHEMICAL CONCENTRATIONS IN THE TOP LEVELS,EXCEPT 
!   FOR AGE
!
!----------------------------------------------------------------------*
!      IF(PRESS(1).LT.8.0) THEN
!        DO IL = 1,MERIDS                                            
!        CH_TEND(IL,21) = VTAUYR     
!       ENDDO
!        RETURN
!     ENDIF
!----------------------------------------------------------------------*
!   RESTRICT THE ACTUAL CHEMICAL MODEL TO STRATOSPHERIC AND 
!   MESOSPHERIC AIR
!----------------------------------------------------------------------*
      DO IL = 1,MERIDS                                               
        DO IZ = 1,NCHEM                                                
          DYY(IL,IZ) = DY(IL,IZ)*RHO(IL)                                    
        ENDDO                  
        H2O(IL) = H2O(IL)*RHO(IL)
        DYZ(IL,1) = OZONE(IL)*RHO(IL)
        DYZ1(IL,1) = OZONE(IL)*RHO(IL)
!        IF(PRESS(IL).LT.10.0.AND.OZONE(IL).GT.1.0E-5) THEN
!          DYZ(IL,1) = 1.0E-5*RHO(IL)
!          DYZ1(IL,1) = DYZ(IL,1)
!        ENDIF                                            
        DO IZ = 2,JIMPL                                                
          DYZ(IL,IZ) = 0.                                                   
          DYZ1(IL,IZ) = 0.                                                  
        ENDDO
      ENDDO                                                   
!----------------------------------------------------------------------*
!   CALL SUBROUTINES FOR TEMPERATURE DEPENDENT RATES, HETEROGENEOUS RATES 
!   AND  FOR PHOTOLYSIS RATES (ONCE PER CALL TO CHEMISTRY)
!   PASS THE PRESSURE IN PA INTO THE OTHER CHEMISTRY S/RS
!----------------------------------------------------------------------*
      CALL REACTION_RATES (RHO,TEMP,VTEMP,PRESS,RATES,MERIDS,IHET,NHET)
      IF(KL.GT.10.AND.KL.LT.28) THEN
        DO IL = 1,MERIDS
!
!  Set H2SO4 as an analytical function of age, peaking at 0.5 ppbv
!
          X1 = 3.0*(ALOG10(PRESS(IL)) - 3.5)**2
          X2 = 3.0/(10.0*DYY(IL,21) + 1.0)*(DYY(IL,21) - 1.0)**2
          H2SO4(IL) = (0.01 + 0.49*EXP(-X1-X2))*1.0E-9*RHO(IL)
        ENDDO
        CALL PSC(TEMP,PRESS,RHO,DYY,H2O,H2SO4,ANAT,AICE,COND,TICE,      &
           WH2SO4,AM,AW,ALIQ,RMEAN,ASAT,RNAT,RICE,MERIDS,NCHEM,MYPE)
        DO IL = 1,MERIDS
          BNAT(IL) = COND(IL,2)/RHO(IL)
          BICE(IL) = COND(IL,3)/RHO(IL)
        ENDDO
        CALL GAMMA_AER(TEMP,PRESS,RHO,DYY,H2O,WH2SO4,AM,AW,RMEAN,GAMMA, &
           RTT,MERIDS,NCHEM,NHET,MYPE)
        CALL HETRATES(TEMP,PRESS,RHO,DYY,H2O,TICE,ANAT,AICE,AEROSOL,    &
           RMEAN,RNAT,RICE,GAMMA,RTT,RATES,IHET,NHET,MERIDS,NCHEM,MYPE)       
      ENDIF
!
      DO NC=1,NCHEM                                                  
        DO IL = 1,MERIDS                                               
          DYY0(IL,NC) = DYY(IL,NC) ! save amounts after natice subtract 
          H2O0(IL) = H2O(IL) 
        ENDDO
      ENDDO                                                
      KLX = KL
      CALL PHOTO_RATES(VTEMP,OZCOL,COZEN,KLX,PHOTO,OZON,COSP,COSPHC,         &
                        MERIDS,IRX,IRX4,FOTO,MYPE)
! ---------------------------------------------------------------------*
!  INCLUDE VARIATIONS IN 11-YR CYCLE AND 27-DAY SOLAR CYCLE AND MULTIPLY 
!  FOTO BY SUN-EARTH DISTANCE TERM (CDX). 
!     IM....NO. OF MONTHS SINCE DEC 1949: USED TO DETERMINE SOLAR F10.7 FLUX.
!     SOLARTIME....TIME AFTER 1 JAN 1960 FOR (ARBITRARY) MIN OF 27 DAY CYCLE
!     ASSUME 12EQUAL MONTHS PER YEAR OF 365.25 DAYS
!  PHOTOLYSIS RATES ARE TAKEN TO BE LINEARLY RELATED TO THE F10.7 FLUX 
!  SOLAR MAX AND SOLAR MIN PHOTO RATES REFER TO APPROX F10.7 FLUX OF 
!  208.1 AND 72.0, CORRESPONDING TO THE YEARS 1991 AND 1996.
! ------------------------------------------------------------------
!
      IF(ISC.EQ.4) THEN
        FACT = 0.0
        SOL = 140.05
      ELSE
        IDAYS = (ITIME(1) - 1960)*365 + IMON(ITIME(2)) + ITIME(3)
        LEAPD = (ITIME(1) - 1960)/4
        JMOD = 4*LEAPD + 1960
        IF(ITIME(2).LT.3.AND.JMOD.EQ.ITIME(1)) LEAPD = LEAPD - 1
        FRAC = FLOAT(ITIME(4)*60 + ITIME(5))/1440.0
        SOLARTIME = FLOAT(IDAYS + LEAPD) + FRAC 
        ICYC = SOLARTIME/27
        SOLAR27 = (SOLARTIME - 27.0*ICYC)/27.0
        FACT = SIN(SOLAR27*2*PI)
        DX1 =  SOLARTIME*12.0/365.25 + 0.5       
        IDX = DX1
        DX1 = DX1 - IDX
        IM = 120 + IDX
        IF(IM.LT.1) IM = 1
        IF(IM.GT.1800) IM = 1800
        SOL = 0.1*(SOLARDATA(IM+1)*DX1 + (1.0 - DX1)*SOLARDATA(IM))
      ENDIF
      DO IR = 1,IRX
        DO IL = 1,MERIDS
          D11 = FOTO(IL,IR) - FOTO(IL,IR+IRX*2)
          D27 = FOTO(IL,IR+IRX) - FOTO(IL,IR+IRX*3)
          F11 = FOTO(IL,IR+IRX*2) + (SOL - 72.0)*D11/136.1
          F27 = FOTO(IL,IR+IRX*3) + (SOL - 72.0)*D27/136.1
          FOTO(IL,IR) = CDX*(F11 + F27*FACT) 
          IF(FOTO(IL,IR).LT.0.0) FOTO(IL,IR) = 0.0
        ENDDO
      ENDDO
! ---------------------------------------------------------------------*
! --- ITX    : TIMESTEP NUMBER                                          
! ---   IT     : ITERATION LOOP COUNTER                                 
! ---                                                                   
! ---   FIRST COMPUTE N + NO + NO2 + NO3 (ANOZ), CL + CLO (CLOX) AND        
! ---   BR + BRO + BRCL (BROX)                                          
! ---                                                                   
! ---------------------------------------------------------------------*
!
      ITX = 0
      ITNO = 6
      ITIME_LEFT = NINT(DT)
      DO 1100
!
      IF(ITIME_LEFT .LE. ISTEP_CHEM) THEN
        DELT = REAL(ITIME_LEFT)
        ITIME_LEFT = 0
      ELSE
        DELT = REAL(ISTEP_CHEM)
        ITIME_LEFT = ITIME_LEFT - ISTEP_CHEM
      ENDIF
      DVDT = 1.0/DELT
      ITX = ITX + 1
      IF (ITX .GE. 2) ITNO = 2
      DO 1000 IT = 1,ITNO
        DO 200 IL = 1,MERIDS
          SUM = DYY(IL,1) + 2.*DYY(IL,2) + DYY(IL,6) + DYY(IL,12) +          &
                DYY(IL,10)
          ANOZ(IL) = DYY(IL,15) - SUM
          IF(ANOZ(IL).LT.0.)  ANOZ(IL) = 0.
          DYZ(IL,5) = ANOZ(IL)
          IF(DYY(IL,1).GT.DYY(IL,15)) DYY(IL,1) = DYY(IL,15)
          SUM = DYY(IL,4) + DYY(IL,5) + DYY(IL,6)
          CLOX(IL) = DYY(IL,16) - SUM
          IF(CLOX(IL).LT.0.)  CLOX(IL) = 0.
          SUM = DYY(IL,12) + DYY(IL,11) + DYY(IL,9)
          BROX(IL) = DYY(IL,17) - SUM
          IF(BROX(IL).LT.0.)  BROX(IL) = 0.
  200   ENDDO
!----------------------------------------------------------------------*
!---
!---         COMPUTE PHOTOCHEMICAL EQUILIBRIUM FOR THE O1D/O3 AND       
!---         O3P/O3 RATIOS
!----------------------------------------------------------------------*
        DO 900 IL=1,MERIDS
          ANUM = (FOTO(IL,2) + FOTO(IL,28)*DYY(IL,20)/DYY(IL,8))
          ADEN =                                                            &
               (RATES(IL,3) + RATES(IL,9)*H2O(IL) +                         &
               RATES(IL,30)*DYY(IL,18) + RATES(IL,66)*DYY(IL,20) +          &
               RATES(IL,70)*DYY(IL,20))
          RATIO2 = ANUM/ADEN
          ANUM = RATES(IL,3)*RATIO2 + FOTO(IL,3) +                          &
                 (2.*FOTO(IL,1)*O2*RHO(IL) + DYZ(IL,5)*FOTO(IL,4) +         &
                  DYZ(IL,4)*FOTO(IL,31) +                                   &
                  DYZ(IL,6)*FOTO(IL,6) + RATES(IL,51)*DYZ(IL,7)**2)/        &
                  DYZ(IL,1)
          ADEN = RATES(IL, 1)           + RATES(IL, 2)*DYZ(IL, 1) +          &
                 RATES(IL, 5)*DYZ(IL,5) + RATES(IL,11)*DYZ(IL, 8) +          &
                 RATES(IL,12)*DYZ(IL,7) + RATES(IL,21)*DYZ(IL,11) +          &
                 RATES(IL,24)*DYY(IL,6) + RATES(IL,44)*DYY(IL, 9)           
          RATIO1 = ANUM/ADEN
          SUM = 1. + RATIO1 + RATIO2
          RAT(IL,1) = 1./SUM
          DYZ(IL,1) = DYY(IL,8)*RAT(IL,1)
          DYZ(IL,2) = RATIO1*DYZ(IL,1)
          DYZ(IL,3) = RATIO2*DYZ(IL,1)
          RAT(IL,2) = RATIO1*RAT(IL,1)
          RAT(IL,3) = RATIO2*RAT(IL,1)
          IF(DYZ(IL,1).LT.1.0E-15*RHO(IL)) THEN
            XYZ = 1.0E-15*RHO(IL) - DYZ(IL,1)
            DYZ(IL,1) = 1.0E-15*RHO(IL)
            DYZ(IL,2) = DYZ(IL,2) + XYZ
            RAT(IL,1) = DYZ(IL,1)/DYY(IL,8)
            RAT(IL,2) = DYZ(IL,2)/DYY(IL,8)
            RAT(IL,3) = 1.0 - RAT(IL,1) - RAT(IL,2)
            IF(RAT(IL,3).LT.0.0) RAT(IL,3) = 0.0
          ENDIF
  900   ENDDO
!----------------------------------------------------------------------*
!--
!---         APPLY PHOTOCHEMICAL EQUILIBRIUM FOR NO3 SUBJECT TO
!---        ENOUGH NOZ BEING PRESENT
!----------------------------------------------------------------------*
        DO 910 IL = 1,MERIDS
          ADEN = RATES(IL,7)*DYZ(IL,5) + FOTO(IL,6) + FOTO(IL,7) +          &
                 RATES(IL,53)*DYZ(IL,4) + RATES(IL,58)*DYZ(IL,5) +          &
                 2.*RATES(IL,59)*DYZ(IL,6) + RATES(IL,62)*DYY(IL,7)            
          IF(ADEN.LT.DVDT.AND.FOTO(IL,6).NE.0.) ADEN = DVDT                 
          IF(ADEN.NE.0.) THEN
            ANUM = (FOTO(IL,8) + RATES(IL,8))*DYY(IL,2) +                   &
              RATES(IL,19)*DYZ(IL, 7)*DYY(IL,1) +                           &
              RATES(IL, 6)*DYZ(IL, 5)*DYZ(IL,1) +                           &
             (RATES(IL,24)*DYZ(IL, 2) + FOTO(IL,12))*DYY(IL,6) +            &
               FOTO(IL,17)*DYY(IL,12)
            DYZ(IL,6) = ANUM/ADEN
            IF(DYZ(IL,6).GT.(0.5*ANOZ(IL))) DYZ(IL,6) = 0.5*ANOZ(IL)
          ENDIF           
!----------------------------------------------------------------------*
!---        APPLY PHOTOCHEMICAL EQUILIBRIUM FOR N SUBJECT TO         
!---       ENOUGH NOZ BEING PRESENT
!---
!----------------------------------------------------------------------*
          ADEN = RATES(IL,67)*DYZ(IL,4) + RATES(IL,68)*O2*RHO(IL) +         &
                 RATES(IL,71)*DYZ(IL,5)
          ANUM = FOTO(IL,31)*DYZ(IL,4)
          DYZ(IL,18) = ANUM/ADEN
          IF(DYZ(IL,18).GT.(0.5*ANOZ(IL))) DYZ(IL,18) = 0.5*ANOZ(IL)
!----------------------------------------------------------------------*
!---
!---    COMPUTE NO/NO2 RATIO ASSUMING PHOTOCHEMICAL EQUILIBRIUM         
!---
! ---   SOLVE SIMULTANEOUSEQUATIONS
! ---      DA1*NO - DA2*NO2 = DD1
! ---    - DB1*NO + DB2*NO2 = DD2
! ---
! ---       DA1*NO......THE RATE OF LOSS TERM OF NO
! ---       DA2*NO2.....THE RATE OF PRODUCTION OF NO FROM NO2
! ---       DB1*NO......THE RATE OF PRODUCTION OF NO2 FROM NO
! ---       DB2*NO2.....THE RATE OF LOSS TERM OF NO2
! ---       DD1.........THE TOTAL PRODUCTION TERM (P) FOR NO LESS
! ---                   CONTRIBUTIONS FROM NO2
! ---       DD2.........THE TOTAL PRODUCTION TERM (P) FOR NO2 LESS
! ---                   CONTRIBUTIONS FROM NO
! ---------------------------------------------------------------------*
          ANOX = ANOZ(IL) - DYZ(IL,6) - DYZ(IL,18)
          DAX = RATES(IL, 4)*DYZ(IL, 1) + RATES(IL,18)*DYZ(IL, 8) +          &
                RATES(IL,22)*DYZ(IL,11) + RATES(IL,43)*DYZ(IL,13) +         &
                RATES(IL,53)*DYZ(IL,6) + RATES(IL,55)*DYZ(IL,16) +          &
                RATES(IL,64)*DYZ(IL,2)          
          DA1 = DAX + RATES(IL,65)*DYZ(IL,7) + RATES(IL,67)*DYZ(IL,18) +    &
                FOTO(IL,31)
          DA2 = FOTO(IL,4) + RATES(IL,5)*DYZ(IL,2) +                        &
                RATES(IL,58)*DYZ(IL,6)
          DD1 = FOTO(IL,7)*DYZ(IL,6) + FOTO(IL,26)*DYZ(IL,17) +             &
                RATES(IL,68)*DYZ(IL,18)*O2*RHO(IL) +                        &
                2.*RATES(IL,66)*DYY(IL,20)*DYZ(IL,3)
          DB1 = DAX + RATES(IL,53)*DYZ(IL,6)
          DB2 = DA2 + RATES(IL, 6)*DYZ(IL, 1) + RATES(IL, 7)*DYZ(IL, 6) +    &
                      RATES(IL,17)*DYZ(IL, 7) + RATES(IL,23)*DYZ(IL,11) +    &
                      RATES(IL,40)*DYZ(IL,13) + RATES(IL,49)*DYZ(IL, 8) +   &
                      RATES(IL,71)*DYZ(IL,18)
          DD2 = RATES(IL, 8)*DYY(IL,2) + FOTO(IL, 5)*DYY(IL,1) +             &
                 FOTO(IL, 6)*DYZ(IL,6) + FOTO(IL, 8)*DYY(IL,2) +             &
               (RATES(IL,50)*DYZ(IL,7) + FOTO(IL,25) + RATES(IL,63))*       &
              DYY(IL,10) + 2.*RATES(IL,59)*DYZ(IL,6)**2 +                   &
              RATES(IL,60)*DYZ(IL,7)*DYZ(IL,17)
          DDNUM = DD2*DA1 + DB1*DD1
          IF(DDNUM.NE.0.) THEN
            RATIO = (DB2*DD1 + DA2*DD2)/DDNUM
          ELSE
            RATIO = 0.
          ENDIF
          DYZ(IL,5) = ANOX/(1. + RATIO)
          DYZ(IL,4) = RATIO*DYZ(IL,5)
  910   ENDDO
!----------------------------------------------------------------------*
!---
!---    COMPUTE THE BR/BRO AND BR/BRCL RATIO FROM THE EQUATIONS         
!---
! ---   SOLVE THREE WAY SIMULTANEOUSEQUATIONS                          
! ---      DA1*BR - DA2*BRO - DA3*BRCL = DD1                            
! ---    - DB1*BR + DB2*BRO            = DD2                            
! ---             - DC2*BRO + DC3*BRCL = DD3                            
! ---
! ---------------------------------------------------------------------*
        DO 930 IL = 1,MERIDS
          DA1 = RATES(IL,39)*DYZ(IL, 1) + RATES(IL,45)*DYZ(IL,8) +           &
                RATES(IL,46)*DYY(IL, 7)
          DA2 = RATES(IL,41)*DYZ(IL,11) + RATES(IL,43)*DYZ(IL,4) +           &
                 FOTO(IL,19)
          DA3 =  FOTO(IL,24)
          DD1 = RATES(IL,    47)*DYZ(IL, 7)*DYY(IL,11) +                     &
             2.*RATES(IL,IHET+7)*DYY(IL,11)*DYY(IL, 9) +                     &
                 FOTO(IL,    17)*DYY(IL,12)            +                     &
                 FOTO(IL,    18)*DYY(IL, 9)
          DB1 = RATES(IL,    39)*DYZ(IL, 1)
          DB2 = DA2 + RATES(IL,40)*DYZ(IL, 5) + RATES(IL,42)*DYZ(IL,8) +     &
                      RATES(IL,48)*DYZ(IL,11)
          DD2 = RATES(IL,44)*DYZ(IL, 2)*DYY(IL,9)
          DC2 = RATES(IL,48)*DYZ(IL,11)
          DC3 = DA3
          DD3 = RATES(IL,IHET+5)*DYY(IL, 9)*DYY(IL,4) +                      &
                RATES(IL,IHET+6)*DYY(IL,11)*DYY(IL,5)
!
          R1 = 0.0
          R2 = 0.0
          DD13 = DD1*DC3 + DD3*DA3
          DDNUM = DD13*DB1 + DD2*DA1*DC3
          IF(DDNUM.GT.0.0) THEN
            VDDNUM = 1.0/DDNUM
            DAC2 = DA2*DC3 + DC2*DA3
            R1 = (DD13*DB2 + DD2*DAC2)*VDDNUM
            DA12 = DA1*DB2 - DA2*DB1
            IF(DA12.LT.0.0) DA12 = 0.0
            DD12 = DA1*DD2 + DB1*DD1
            R2 = (DD12*DC2 + DD3*DA12)*VDDNUM
          ENDIF
          IF (FOTO(IL,24).GT.0.0) THEN    
            DYZ(IL,13) = BROX(IL)/(1. + R1 + R2)
            DYZ(IL,12) = R1*DYZ(IL,13)
            DYZ(IL,15) = R2*DYZ(IL,13)
          ELSE
            DYZ(IL,13) = 0.0
            DYZ(IL,12) = 0.0
            DYZ(IL,15) = BROX(IL)
          ENDIF
  930   ENDDO
!----------------------------------------------------------------------*
!---
!---    COMPUTE CL/CLO RATIO ASSUMING PHOTOCHEMICAL EQUILIBRIUM         
!---
! ---   SOLVE SIMULTANEOUS EQUATIONS
! ---      DA1*CL - DA2*CLO = DD1
! ---    - DB1*CL + DB2*CLO = DD2
! ---
! ---------------------------------------------------------------------*
        DO 940 IL = 1,MERIDS
          DA1 = RATES(IL,20)*DYZ(IL, 1) + RATES(IL,25)*DYY(IL,18) +         &
                RATES(IL,27)*DYZ(IL, 8) + RATES(IL,33)*DYY(IL, 7)      
          DA2 = RATES(IL,21)*DYZ(IL, 2) + RATES(IL,22)*DYZ(IL, 4) +         &
                RATES(IL,29)*DYZ(IL, 7) + RATES(IL,41)*DYZ(IL,13)      
          DD1 = RATES(IL,26)*DYZ(IL, 7)*DYY(IL,4) + FOTO(IL,15)*DYY(IL,5) +  &
              2.*FOTO(IL,16)*DYZ(IL,14) +  FOTO(IL,12)*DYY(IL,6) +           &
                 FOTO(IL,24)*DYZ(IL,15) +  FOTO(IL,33)*DYY(IL,4)
          DB1 = RATES(IL,20)*DYZ(IL, 1)
          DB2 = RATES(IL,21)*DYZ(IL, 2) + RATES(IL,22)*DYZ(IL, 4)    +       &
                RATES(IL,23)*DYZ(IL, 5) + RATES(IL,28)*DYZ(IL, 8)    +       &
               (RATES(IL,29) + RATES(IL,69))*DYZ(IL, 7) +                   &
                RATES(IL,37)*DYZ(IL,11)*2. +                                &
               (RATES(IL,41) + RATES(IL,48))*DYZ(IL,13)
          DD2 = RATES(IL,24)*DYZ(IL, 2)*DYY(IL,6) +                          &
             2.*RATES(IL,38)*DYZ(IL,14)
          DDNUM = DD2*DA1 + DB1*DD1
          IF(DDNUM.NE.0.) THEN
            RATIO = (DB2*DD1 + DA2*DD2)/DDNUM
          ELSE
            RATIO = 0.
          ENDIF
          CLX = CLOX(IL) - DYZ(IL,15)
          IF(CLX.LT.0.0) CLX = 0.0
          IF(FOTO(IL,15).GT.0.0)  THEN
            TAU = RATES(IL,37)/(FOTO(IL,16) + RATES(IL,38))                
            RAT1 = TAU*DYZ(IL,11)
            DYZ(IL,11) = CLX/(1. + RATIO + 2.*RAT1)
            RAT1 = TAU*DYZ(IL,11)
            DYZ(IL,11) = CLX/(1. + RATIO + 2.*RAT1)
            DYZ(IL,10) = DYZ(IL,11)*RATIO
            DYZ(IL,14) = DYZ(IL,11)*RAT1
          ELSE
            DYZ(IL,10) = 0.
            DYZ(IL,11) = 0.
            DYZ(IL,14) = CLX*0.5
          ENDIF
  940   ENDDO
!----------------------------------------------------------------------*
!---
!---    COMPUTE HO2 AND OH ASSUMING PHOTOCHEMICAL EQUILIBRIUM           
!---
! ---   SOLVE SIMULTANEOUSEQUATIONS
! ---      DA1*OH - DA2*HO2 = DD1
! ---    - DB1*OH + DB2*HO2 = DD2
! ---
! ---------------------------------------------------------------------*
        DO 950 IL = 1,MERIDS
          BX1 = RATES(IL,17)*DYZ(IL, 5) +                                   &
                RATES(IL,19)*DYY(IL, 1) + RATES(IL,26)*DYY(IL, 4) +         &
                RATES(IL,31)*DYY(IL,18) + RATES(IL,50)*DYY(IL,10) +         &
                RATES(IL,57)*DYY(IL,13) + RATES(IL,60)*DYZ(IL,17) +         &
                RATES(IL,65)*DYZ(IL,4) + RATES(IL,69)*DYZ(IL,11)
          BX2 = RATES(IL,27)*DYZ(IL,10) + RATES(IL,28)*DYZ(IL,11) +         &
                RATES(IL,42)*DYZ(IL,13) + RATES(IL,45)*DYZ(IL,12) +         &
                RATES(IL,49)*DYZ(IL, 5) + RATES(IL,56)*DYZ(IL,16)
          DB1 = RATES(IL,10)*DYZ(IL, 1) + RATES(IL,12)*DYZ(IL, 2) +         &
                RATES(IL,29)*DYZ(IL,11) + RATES(IL,32)*DYY(IL, 7) +         &
                RATES(IL,36)*DYY(IL,3) +  RATES(IL,47)*DYY(IL,11) +         &
                RATES(IL,52)*DYY(IL,14) + RATES(IL,54)*H2*RHO(IL)
          DA1 = DB1 + BX1 + RATES(IL,16)*DYZ(IL, 8) +                       &
                 (RATES(IL,35) + RATES(IL,51))*DYZ(IL,7)*2.
          DA2 = RATES(IL,11)*DYZ(IL, 2) + RATES(IL,15)*DYZ(IL,1) +          &
                RATES(IL,18)*DYZ(IL, 4)
          DB2 = RATES(IL,16)*DYZ(IL, 7) + RATES(IL,34)*DYZ(IL,8)*2. +       &
                DA2 + BX2
          DD1 = DYZ(IL,3)*(2.*RATES(IL, 9)*H2O(IL) +                        &
                              RATES(IL,30)*DYY(IL,18)) +                    &
                RATES(IL,44)*DYZ(IL,2)*DYY(IL,9) +                          &
                2.*RATES(IL,61)*DYZ(IL,16)**2 +                             &
                 FOTO(IL, 5)*DYY(IL, 1) +  FOTO(IL,15)*DYY(IL, 5) +         &
              2.*FOTO(IL,11)*DYY(IL, 3) +  FOTO(IL,18)*DYY(IL, 9) +         &
                FOTO(IL,26)*DYZ(IL,17) +                                    &
! FOTO(IL,27) IS TAKEN AS A PROXY FOR CH3O2 + HV --> H2CO + OH      
                FOTO(IL,27)*(DYZ(IL,16)+DYY(IL,13)) +                       &
                FOTO(IL,32)*H2O(IL)
          DD2 = (RATES(IL,33)*DYZ(IL,10) +                                  &
                RATES(IL,46)*DYZ(IL,12) + FOTO(IL,13)*2.0)*DYY(IL,7) +      &
                RATES(IL,55)*DYZ(IL,4)*DYZ(IL,16) +                         &
                RATES(IL,62)*DYZ(IL,6)*DYY(IL,7) +                          &
                RATES(IL,63)*DYY(IL,10) + FOTO(IL,25)*DYY(IL,10) +          &
                FOTO(IL,27)*DYY(IL,13) +                                    &
                FOTO(IL,32)*H2O(IL)  +  FOTO(IL,33)*DYY(IL,4) +             &
                FOTO(IL,32)*DYY(IL,18)*1.53/1.83
          DDNUM = DB2*DD1 + DA2*DD2
          IF(DD1.NE.0.0.AND.DD2.NE.0.0) THEN
            RATIO2 = (DD2*DA1 + DB1*DD1)/DDNUM
            AX = 4.*(RATES(IL,16)*RATIO2 + RATES(IL,34)*RATIO2**2 +        &
                  RATES(IL,35) + RATES(IL,51))
            BX = BX1 + RATIO2*BX2
            CX = DD1 + DD2
            IF(CX.LT.0.) CX = 0.0
            IF(BX.LT.0.) BX = 0.0
            IF(BX.EQ.0.) THEN
              DYZ(IL,7) = SQRT(2.*CX/AX)
            ELSE
              DYZ(IL,7) = (-BX + SQRT(BX**2 + 2.*AX*CX))/AX
            ENDIF
            DYZ(IL,8) = RATIO2*DYZ(IL,7)
          ELSE
            DYZ(IL,7) = 0.0
            DYZ(IL,8) = 0.0
          ENDIF
          RATIO1 = RATES(IL,12)*DYZ(IL,2)/(RATES(IL,14) +                   &
                   RATES(IL,13)*DYZ(IL,1))
          DYZ(IL,9) = RATIO1*DYZ(IL,7)
  950   ENDDO
!----------------------------------------------------------------------*
!---
!---    COMPUTE CH3O2 AND HONO ASSUMING PHOTOCHEMICAL EQUILIBRIUM           
!---
!----------------------------------------------------------------------*
        DO 960 IL = 1,MERIDS
          AX = 4.*RATES(IL,61)
!
! FOTO(IL,27) IS TAKEN AS A PROXY FOR CH3O2 + HV --> H2CO + OH
!
          BX = RATES(IL,55)*DYZ(IL,4) + RATES(IL,56)*DYZ(IL,8) +            &
               FOTO(IL,27)
          CX = (RATES(IL,25)*DYZ(IL,10) + RATES(IL,30)*DYZ(IL,3) +          &
                RATES(IL,31)*DYZ(IL,7) + FOTO(IL,32)*1.53/1.83)*DYY(IL,18) + &
                RATES(IL,57)*DYZ(IL,7)*DYY(IL,13)
          IF(CX.LT.0.) CX = 0.0
          IF(BX.LT.0.) BX = 0.0
          IF(BX.EQ.0.) THEN
            DYZ(IL,16) = SQRT(2.*CX/AX)
          ELSE
            DYZ(IL,16) = (-BX + SQRT(BX**2 + 2.*AX*CX))/AX
          ENDIF
  960   ENDDO!CONTINUE
        DO 970 IL = 1,MERIDS
          DDNUM = RATES(IL,60)*DYZ(IL,7) + FOTO(IL,26)
          IF(DDNUM.NE.0.0)  THEN
            DYZ(IL,17) = RATES(IL,65)*DYZ(IL,7)*DYZ(IL,4)/DDNUM
          ELSE
            DYZ(IL,17) = RATES(IL,65)*DYZ(IL,4)/RATES(IL,60)
          ENDIF
  970   ENDDO!CONTINUE
        IF(IT.EQ.1.AND.ITX.EQ.1) THEN
          DO IC = 1,JIMPL
            DO IL = 1,MERIDS
              DYZ1(IL,IC) = DYZ(IL,IC)
            ENDDO
          ENDDO
        ELSE
          DO IC = 1,JIMPL
            DO IL = 1,MERIDS
              XYZ = (DYZ(IL,IC) + DYZ1(IL,IC))*0.5
              DYZ(IL,IC) = XYZ
              DYZ1(IL,IC) = XYZ
            ENDDO
          ENDDO
        ENDIF
!----------------------------------------------------------------------*
!   INTEGRATE MODEL FOR SINGLE TIMESTEP OF LENGTH DT
! --- THE CHEMISTRY SCHEME CALCULATES DESTRUCTION AND
! --- PRODUCTION TERMS READY FOR AN INTEGRATION OF THE FORM
! ---
! ---                 Y = P/Q + (Y - P/Q)*EXP(-Q.DT)
! ---
! ---EXCEPT WHERE Q.DT IS SMALL, WHERE 
! ---   
! ---                 Y = (Y/DT + P)/(Q + 1/DT)
! ---
! --- (FIRST ORDER EULER BACKWARD SCHEME) IS USED.
! ---
! ---   DQT : TOTAL DESTRUCTION RATE (FOR EACH EXPLICIT SPECIES)        
! ---   DPT : TOTAL PRODUCTION RATE (FOR EACH EXPLICIT SPECIES)         
! ---
! ---------------------------------------------------------------------*
        DO 990 IL = 1,MERIDS
          JLX = 1 + (PI*0.5 + ALAT(IL))*89.0/PI
!
!   HNO3
!
          DPT(IL,1) = RATES(IL,17)*DYZ(IL,5)*DYZ(IL,7) +                    &
                      RATES(IL,62)*DYZ(IL,6)*DYY(IL,7) +                    &
                      RATES(IL,IHET+1)*DYY(IL, 2)*DYY(IL, 4)    +           &
                      RATES(IL,IHET+2)*H2O(IL)*DYY(IL, 2)*2. +              &
                      RATES(IL,IHET+3)*H2O(IL)*DYY(IL, 6)    +              &
                      RATES(IL,IHET+4)*DYY(IL, 6)*DYY(IL, 4)
!
!  The HNO3 production term from heterogeneous reaction BrONO2 + H2O
!  is sometimes rather fast leading to too much HNO3 extrapolated to the 
!  full timestep. Hence this term is neglected since it should be small in 
!  any case.      
!     *            RATES(IL,IHET+8)*DYY(IL,12)*H2O(IL)  
!
          DQT(IL,1) = RATES(IL,19)*DYZ(IL,7) + FOTO(IL,5)            
!       
!   N2O5
!
          DPT(IL,2) = RATES(IL,7)*DYZ(IL,5)*DYZ(IL,6)
          DQT(IL,2) = RATES(IL,8) + FOTO(IL,8) +                            &
                      RATES(IL,IHET+1)*DYY(IL, 4) +                         &
                      RATES(IL,IHET+2)*H2O(IL)
!
!  H2O2
!
          DPT(IL,3) = RATES(IL,34)*DYZ(IL,8)**2 +                           &
                      RATES(IL,35)*DYZ(IL,7)**2
          DQT(IL,3) = FOTO(IL,11) + RATES(IL,36)*DYZ(IL,7)
!
!  HCl
!
          DPT(IL,4) = DYZ(IL,10)*(RATES(IL,27)*DYZ(IL, 8)  +                &
                                  RATES(IL,33)*DYY(IL, 7)  +                &
                                  RATES(IL,25)*DYY(IL,18)) +                &
                      RATES(IL,69)*DYZ(IL,11)*DYZ(IL,7)
          DQT(IL,4) = RATES(IL,    26)*DYZ(IL, 7)    +                      &
                      RATES(IL,  IHET)*DYY(IL, 5)    +                      &
                      RATES(IL,IHET+1)*DYY(IL, 2)    +                      &
                      RATES(IL,IHET+4)*DYY(IL, 6)    +                      &
                      RATES(IL,IHET+5)*DYY(IL, 9)    +  FOTO(IL,33)
!
!  HOCl
!
          DPT(IL,5) = RATES(IL,    28)*DYZ(IL,11)*DYZ(IL, 8) +              &
                      RATES(IL,IHET+3)*DYY(IL, 6)*H2O(IL)
          DQT(IL,5) = FOTO(IL,15)                 +                         &
                      RATES(IL,  IHET)*DYY(IL, 4) +                         &
                      RATES(IL,IHET+6)*DYY(IL,11)
!
!  ClONO2
!
          DPT(IL,6) = RATES(IL,23)*DYZ(IL,11)*DYZ(IL,5)
          DQT(IL,6) = RATES(IL,24)*DYZ(IL, 2) + FOTO(IL,12) +               &
                      RATES(IL,IHET+3)*H2O(IL) +                            &
                      RATES(IL,IHET+4)*DYY(IL, 4)
!
!  H2CO
!
          DPT(IL,7) = DYZ(IL,16)*(RATES(IL,55)*DYZ(IL,4) +                  &
                    2.*RATES(IL,61)*DYZ(IL,16)) + FOTO(IL,27)*DYY(IL,13)
          DQT(IL,7) = FOTO(IL,13) + FOTO(IL,14) +                           &
                      RATES(IL,32)*DYZ(IL, 7) + RATES(IL,33)*DYZ(IL,10) +   &
                      RATES(IL,46)*DYZ(IL,12) + RATES(IL,62)*DYZ(IL,6)          
!
! Oy
!
          DPT(IL,8) = RATES(IL,51)*DYZ(IL,7)**2 +                           &
                   2.*FOTO(IL,1)*O2*RHO(IL) + FOTO(IL, 4)*DYZ(IL, 5) +      &
                      FOTO(IL,19)*DYZ(IL,13) + FOTO(IL,28)*DYY(IL,20)    
          DQT1 = RAT(IL,1)*                                                 &
                (RATES(IL, 4)*DYZ(IL, 4) +                                  &
                RATES(IL,10)*DYZ(IL, 7) + RATES(IL,13)*DYZ(IL, 9) +         &
                RATES(IL,15)*DYZ(IL, 8) + RATES(IL,20)*DYZ(IL,10) +         &
                RATES(IL,39)*DYZ(IL,12))
!     &      RATES(IL,15)*DYZ(IL, 8) + RATES(IL,20)*DYZ(IL,10))
          DQT2 = RAT(IL,2)*                                                 &
            (2.*RATES(IL, 2)*DYZ(IL, 1) + RATES(IL, 5)*DYZ(IL, 5) +         &
                RATES(IL,11)*DYZ(IL, 8) + RATES(IL,12)*DYZ(IL, 7) +         &
                RATES(IL,21)*DYZ(IL,11) + RATES(IL,24)*DYY(IL, 6) +         &
                RATES(IL,44)*DYY(IL,9)  + RATES(IL,64)*DYZ(IL, 4))
          DQT3 = RAT(IL,3)*                                                 &
               (RATES(IL, 9)*H2O(IL) +                                      &
                RATES(IL,30)*DYY(IL,18))
!      DQT4 = RATES(IL,39)*OZONE(IL)*ISTEP_CHEM
!      IF(DQT4.GT.1.0) DQT4 = 1.0 
!      DQT(IL,8) = DQT1 + DQT2 + DQT3 + DQT4*DYZ(IL,12)/(ISTEP_CHEM*OZONE(IL))
          DQT(IL,8) = DQT1 + DQT2 + DQT3 
!
!  HOBr
!
          DPT(IL,9) = RATES(IL,    42)*DYZ(IL,13)*DYZ(IL, 8) +              &
                      RATES(IL,IHET+8)*H2O(IL)*DYY(IL,12)
          DQT(IL,9) = FOTO(IL,18) + RATES(IL,44)*DYZ(IL,2) +                &
                      RATES(IL,IHET+5)*DYY(IL, 4)          +                &
                      RATES(IL,IHET+7)*DYY(IL,11)
!  
!  HNO4
!
          DPT(IL,10) = RATES(IL,49)*DYZ(IL,8)*DYZ(IL,5)
          DQT(IL,10) = FOTO(IL,25) + RATES(IL,50)*DYZ(IL,7) +               &
                       RATES(IL,63)
!
!  HBr
!
          DPT(IL,11) = RATES(IL,45)*DYZ(IL,12)*DYZ(IL,8) +                  &
                       RATES(IL,46)*DYZ(IL,12)*DYY(IL,7)
          DQT(IL,11) = RATES(IL,    47)*DYZ(IL, 7) +                        &
                       RATES(IL,IHET+6)*DYY(IL, 5) +                        &
                       RATES(IL,IHET+7)*DYY(IL, 9)
!
!  BrONO2
!
          DPT(IL,12) = RATES(IL,40)*DYZ(IL,13)*DYZ(IL,5)
          DQT(IL,12) = FOTO(IL,17) + RATES(IL,IHET+8)*H2O(IL)
!
!  CH3OOH
!
          DPT(IL,13) = RATES(IL,56)*DYZ(IL,16)*DYZ(IL,8)
          DQT(IL,13) = FOTO(IL,27) + RATES(IL,57)*DYZ(IL,7)
!
! CO
!
          DPT(IL,14) = DYY(IL,7)*                                           &
                  (RATES(IL,32)*DYZ(IL,7) + RATES(IL,33)*DYZ(IL,10) +       &
                   RATES(IL,46)*DYZ(IL,12) + RATES(IL,62)*DYZ(IL,6) +       &
                   FOTO(IL,13) + FOTO(IL,14))
          DQT(IL,14) = RATES(IL,52)*DYZ(IL,7)
!
!  NOY + NAT: Normal Stratospheric terms, relax towards ANOY in the 
!  troposphere
!
          AGESQ = 1.0E-2
          IF(DAGESQ(IL).LT.AGESQ) THEN
            DPT(IL,15) = ANOY(JLX,KL)*RHO(IL)*VTAU10
            DQT(IL,15) = VTAU10
          ELSE
            DPT(IL,15) = 2.*RATES(IL,66)*DYZ(IL,3)*DYY(IL,20)
            DQT(IL,15) = 2.*DYZ(IL,18)*                                     &
                   (RATES(IL,67)*DYZ(IL,4)+RATES(IL,71)*DYZ(IL,5))/DYY(IL,15)
          ENDIF
!
!  Cly: solved separately after call to chemistry; set to zero here
!
          DPT(IL,16) = 0.0
          DQT(IL,16) = 0.0
!
!  Bry: solved separately after call to chemistry; set to zero here
!
          DPT(IL,17) = 0.0
          DQT(IL,17) = 0.0
!
!  CH4: Production rate assumes relaxation to a uniform tropospheric value 
!  with a 10 day time scale
!  DAGESQ < 0.01 is used to denote the troposphere. 
!
          IF(DAGESQ(IL).LT.AGESQ) THEN
            DPT(IL,18) = CH4*VTAU10*RHO(IL)
            DQT(IL,18) = VTAU10 
          ELSE
            DPT(IL,18) = 0.0
!  
! Mesospheric loss of CH4 by lyman alpha photolysis is taken to be the H2O 
! photolysis rate weighted by the ratio of the cross sections
! 
            DQT(IL,18) = RATES(IL,25)*DYZ(IL,10) + RATES(IL,30)*DYZ(IL,3) +   &
                         RATES(IL,31)*DYZ(IL,7) + FOTO(IL,32)*1.53/1.83
          ENDIF
!
!  h2ostrat: not used
!
          DPT(IL,19) = 0.0
          DQT(IL,19) = 0.0
!
!  N2O: Production rate assumes relaxation to a uniform tropospheric value 
!  with a 10 day time scale
!  DAGESQ < 0.01 is used to denote the troposphere. 
!
          IF(DAGESQ(IL).LT.AGESQ) THEN
            DPT(IL,20) = AN2O*VTAU10*RHO(IL)
            DQT(IL,20) = VTAU10
          ELSE
            DPT(IL,20) = RATES(IL,71)*DYZ(IL,5)*DYZ(IL,18)
            DQT(IL,20) = (RATES(IL,66)+RATES(IL,70))*DYZ(IL,3) +           &
                          FOTO(IL,28)
          ENDIF
!
!  AGE OF AIR: Increase at the rate of 1 per second outside the troposphere.
!  Results expressed in years. Relax towards zero with a 10 
!  day timescale in the troposphere, denoted by DAGESQ less than 0.01 
!
          IF(DAGESQ(IL).LT.AGESQ) THEN
            DPT(IL,21) = 0.0
            DQT(IL,21) = VTAU10
          ELSE
            DPT(IL,21) = VTAUYR*RHO(IL)
            DQT(IL,21) = 0.0
          ENDIF
!
!  H2O + ICE
!
          DP1(IL) = DYZ(IL,7)*( RATES(IL,16)*DYZ(IL,8) +                    &
             RATES(IL,19)*DYY(IL,1) + RATES(IL,26)*DYY(IL,4) +              &
             RATES(IL,31)*DYY(IL,18) + RATES(IL,32)*DYY(IL,7) +             &
             RATES(IL,36)*DYY(IL,3) + RATES(IL,50)*DYY(IL,10) +             &
             RATES(IL,47)*DYY(IL,11) +                                      &
             RATES(IL,51)*DYZ(IL,7) + RATES(IL,54)*H2*RHO(IL) +             &
             RATES(IL,57)*DYY(IL,13) + RATES(IL,60)*DYZ(IL,17) )
          DQ1(IL) = FOTO(IL,32) + RATES(IL,9)*DYZ(IL,3)
  990   ENDDO
! ---------------------------------------------------------------------*
! ---
! ---  INTEGRATE NCHEM SPECIES EXPLICITLY
! ---
! ---------------------------------------------------------------------*
        DO 980 NC = 1,NCHEM
          DO  IL = 1,MERIDS
            T0 = DQT(IL,NC)/DVDT
            IF(T0 .LT. 1.0E-6) THEN 
              DYY1(IL,NC) = (DYY(IL,NC)*DVDT + DPT(IL,NC))/(DQT(IL,NC)+DVDT)
! The following line changes answers but should be equivalent to the above..
!              DYY1(IL,NC) = (DYY(IL,NC) + DPT(IL,NC)/DVDT)/(T0+1.0)
            ELSE
              minusT0 = -1.0*T0
              T1 = DPT(IL,NC)/DQT(IL,NC)
!              DYY1(IL,NC) = T1 + (DYY(IL,NC) - T1)*EXP(-1.0*T0)
              DYY1(IL,NC) = T1 + (DYY(IL,NC) - T1)*EXP(minusT0)
            ENDIF
          ENDDO
  980   ENDDO!CONTINUE 
        DO  IL = 1,MERIDS
          H2O1(IL) = (H2O(IL)*DVDT + DP1(IL))/(DQ1(IL)+DVDT)
!
! SPECIAL CODE TO ENSURE THAT THE BROMINE SPECIES ARE UNDER CONTROL DURING
! THE RAPIDLY VARYING PERIOD AT DAWN AND DUSK
!
          BRX = DYY1(IL,9) + DYY1(IL,11) + DYY1(IL,12) + DYZ(IL,12) +       &
                DYZ(IL,13) + DYZ(IL,15)
          IF(DYY1(IL,17).EQ.0.0) THEN
            RATIO = 0.0
          ELSE
            RATIO = BRX/DYY1(IL,17)
          ENDIF
          IF(RATIO.GT.1.0) THEN
            DYY1(IL,9) = DYY1(IL,9)/RATIO
            DYY1(IL,11) = DYY1(IL,11)/RATIO
            DYY1(IL,12) = DYY1(IL,12)/RATIO
            DYZ(IL,12) = DYZ(IL,12)/RATIO
            DYZ(IL,13) = DYZ(IL,13)/RATIO
            DYZ(IL,15) = DYZ(IL,15)/RATIO
          ENDIF
!
!  SPECIAL CODE TO ENSURE THAT THE HETEROGENEOUS SCHEME DOES NOT PRODUCE
! EXCESSIVE HNO3. IF HNO3 EXCEEDS NOy - 2.*N2O5 SET HNO3 TO  NOy - 2.*N2O5
!
          DIFF = DYY1(IL,15) - 2.0*DYY1(IL,2)
          IF (DYY1(IL,1).GT.DIFF) DYY1(IL,1) = DIFF
        ENDDO
 1000 ENDDO
!
!  COMPUTE INCREMENT TO O3_CHEM DIAGNOSTIC TRACER
!
      DO IL = 1,MERIDS
        O3_PROD(IL) = O3_PROD(IL) +                                     &
             (DPT(IL,8) - DQT(IL,8)*DYY(IL,8))*RAT(IL,1)*DELT/RHO(IL)
      ENDDO
      IF(ITIME_LEFT .GT. 0) THEN
        DO NC = 1,NCHEM
          DO IL = 1,MERIDS
            DYY(IL,NC) = DYY1(IL,NC)
            H2O(IL) = H2O1(IL)
          ENDDO
        ENDDO
        CYCLE
      ELSE
        EXIT
      ENDIF
 1100 ENDDO
! ---------------------------------------------------------------------*
! ---
! --- COMPUTE CHEMICAL INCREMENTS FOR THIS SUPER-TIMESTEP   
! ---
! ---------------------------------------------------------------------*
!
!  GET CHEM TENDENCY AND CONVERT TO MIXING RATIOS
!
      DO IL = 1,MERIDS
        VRHO = 1.0/(RHO(IL)*DT)
        DO NC = 1,NCHEM
          CH_TEND(IL,NC) = CH_TEND(IL,NC) +                             &
                           (DYY1(IL,NC) - DYY0(IL,NC))*VRHO
       ENDDO
          H2O_TEND(IL) = H2O_TEND(IL) + (H2O1(IL) - H2O0(IL))*VRHO
      ENDDO
!
! SET OZONE CONCENTRATION
!
      DO IL = 1,MERIDS
        OZONE(IL) = DYZ(IL,1)/RHO(IL)
      ENDDO
!
#ifdef _ALLOC
      deallocate(RATES)
      deallocate(FOTO)
      deallocate(DYY)
      deallocate(DYZ)
      deallocate(DQT)
      deallocate(DPT)
      deallocate(RAT)
      deallocate(DYZ1)
      deallocate(DYY0)
      deallocate(DYY1)
      deallocate(COND)
      deallocate(GAMMA)
#endif
      RETURN
END SUBROUTINE CHEMISTRY
!         
SUBROUTINE DCLY_DT(age,dfdage,tropc,tracer1,cly,bry,dclydt,dbrydt, &
           agefact2,merids,lats,levels,itime)       
!----------------------------------------------------------------------*
!      age:  Age of air since stratospheric entry in years.                
!
!      The order of dfdage array elements is:
!
!      1. CFC11 as a fraction of its tropospheric concentration.        
!      2.-7.  CFC12, CFC113, CCL4, CH3CL, CH3CCL3, HCFC22 
!      8  Bry as a fraction of its tropospheric concentration
!
!      TROPC contains the tropospheric concentrations for 1950.0-2100.0 
!      at annual intervals. The species are ordered as dfdage elements 1-7 
!      and Total organic Cl and Br.
!----------------------------------------------------------------------*
      IMPLICIT NONE
      
      integer, intent(in) :: MERIDS,LATS,LEVELS,itime(6) 
      REAL, intent(in) :: age(merids,lats,levels),dfdage(90,48,8),TROPC(151,9), &
            cly(merids,lats,levels),bry(merids,lats,levels),         &
            tracer1(merids,lats,levels)
      REAL, intent(inout) :: dclydt(merids,lats,levels), &
            dbrydt(merids,lats,levels)
      REAL, intent(in) :: AGEFACT2
            
      
      INTEGER MDAY(12)
      INTEGER IANN,IMON,IDAY,IMON1,IT1,IT2,IC,KL,JL, &
              IL,IY,IY1,IM,MX
      REAL  CLWEIGHT(7),                    &
            dfdtau(merids,lats,8),          &
            cfc(merids,lats,9),cfct(9,2)

      REAL  YSTART,TFACT,AGEFACT1,DX1,DX2,TIME,DT1,DT2,SUM1,   &
            SUM2,FACTOR,YY,Y1,Y2,clytot,brytot
      DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA YSTART/1950.0/
      DATA CLWEIGHT/ 3.0, 2.0, 3.0, 4.0, 1.0, 3.0, 1.0 / 
!
!      Interpolate between months for TAB array
!
      iann = itime(1)
      imon = itime(2)
      iday = itime(3)
      tfact = 1.0/(365.25*86400.0)
!
! agefact1 and agefact2 are fudge factors, to correct for low model 
! atmospheric age in the computation of the halogen rates of change.
!
      agefact1 = 1.25
      DX1 = IDAY/MDAY(IMON)
      IF(IMON.EQ.2.AND.4*(IANN/4).EQ.IANN) DX1 = IDAY/29.0
      DX2 = 1.0 - DX1
      IMON1 = IMON + 1
      IF(IMON1.EQ.13) IMON1 = 1
!
!  Compute multiplying factor for missing CFCs, and include factor for 
!  conversion of rates to a per second rate. 
!
      TIME = REAL(IANN) + (REAL(imon)-1.0)/12.0
      IT1 = INT(TIME - YSTART) + 1
      IF(IT1 .LT. 1) IT1 = 1
      IF(IT1 .GT. 150) IT1 = 150
      IT2 = IT1 + 1
      DT1 = TIME - REAL(IT1) - (YSTART - 1.0)
      DT2 = 1.0 - DT1
      sum1 = 0.0
      sum2 = 0.0
      do ic = 1,7
        sum1 = sum1 + clweight(ic)*tropc(it1,ic)
        sum2 = sum2 + clweight(ic)*tropc(it2,ic)
      enddo
      factor = (dt2*tropc(it1,8) + dt1*tropc(it2,8))*tfact/             &
                (sum1*dt2 + sum2*dt1)
      do 500 kl = 1,levels
!----------------------------------------------------------------------*
!  determine age and dfdage at air parcel positions using tracer
!----------------------------------------------------------------------*
        do jl = 1,lats
          do il = 1,merids
            yy = tracer1(il,jl,kl)
            IF(YY.LT.1.)  YY =  1.01
            IF(YY.GT.89.99) YY = 89.99
            IY = YY
            Y1 = YY - IY
            Y2 = 1. - Y1
            IY1 = IY + 1
            do ic = 1,8
            dfdtau(il,jl,ic) = y2*dfdage(iy,kl,ic) + y1*dfdage(iy1,kl,ic)
            enddo
          enddo
        enddo
!----------------------------------------------------------------------*
!  compute CFCs at time t - age
!----------------------------------------------------------------------*
!      DO 100 JL = 1,lats
!      DO 100 IL = 1,merids
        DO JL = 1,lats
          DO IL = 1,merids
            DO IM = 1,2
              MX = IM + IMON - 1
              IF(MX.EQ.13) MX = 1
              TIME = REAL(IANN) + (REAL(MX)-1.0)/12.0 - age(IL,JL,kl)*agefact1  
!
!      Interpolate the tropospheric gas concs in time
!
              IT1 = INT(TIME - YSTART) + 1
              IF(IT1 .LT. 1) IT1 = 1
              IF(IT1 .GT. 150) IT1 = 150
              IT2 = IT1 + 1
              DT1 = TIME - REAL(IT1) - (YSTART - 1.0)
              DT2 = 1.0 - DT1
              DO IC=1,9
                CFCT(IC,IM) = TROPC(IT1,IC)*DT2 + TROPC(IT2,IC)*DT1                
              ENDDO! IC
            ENDDO! IM
            DO IC = 1,9
              CFC(IL,JL,IC) =  CFCT(IC,1)*DX2 +  CFCT(IC,2)*DX1  
            ENDDO !IC
          ENDDO
        ENDDO
! 100  CONTINUE
!----------------------------------------------------------------------*
!  Finally compute dclydt and dbrydt
!----------------------------------------------------------------------*
        do jl = 1,lats
          do il = 1,merids
            dclydt(il,jl,kl) = 0.0
            do ic = 1,7
              dclydt(il,jl,kl) = dclydt(il,jl,kl) +             &
                  1.0e-12*factor*dfdtau(il,jl,ic)*clweight(ic)* &
                  cfc(il,jl,ic)*agefact2
            enddo
            clytot = cfc(il,jl,8)*1.0e-12
            if(cly(il,jl,kl).ge.clytot) dclydt(il,jl,kl) = 0.0
            dbrydt(il,jl,kl) = 1.0e-12*tfact*dfdtau(il,jl,8)*cfc(il,jl,9)
            brytot = 1.0e-12*cfc(il,jl,9)
            if(bry(il,jl,kl).ge.brytot) dbrydt(il,jl,kl) = 0.0
          enddo
        enddo
 500  ENDDO
      RETURN
END SUBROUTINE DCLY_DT
!
SUBROUTINE REACTION_RATES (RHO,TEMP,VTEMP,PRESS,RATES,            &
                           MERIDS,IHET,NHET)
!-------------------------------------------------------------------------
!
!   THIS SUBROUTINE DETERMINES THE HOMOGENEOUS GAS PHASE REACTION RATES
!
!-------------------------------------------------------------------------
IMPLICIT NONE
INTEGER  , intent(in) :: MERIDS,IHET,NHET
REAL, intent(in), dimension(MERIDS) ::  RHO,TEMP,VTEMP,PRESS
REAL, intent(inout) :: RATES(MERIDS,IHET+NHET-1)
!
INTEGER  IL,IC
REAL TEMP3(MERIDS)
!
REAL  O2,R0,R1,RX,RX1,RX2
!
  DATA O2/0.2095/      
!
  DO IL = 1,MERIDS
    TEMP3(IL) = TEMP(IL)/300.0
  ENDDO   
!
!  RATES TAKEN FROM JPL 2006 
!
  DO IL = 1,MERIDS
!   REACTION 1 O+O2+M --> O3+M
      RATES(IL,1) = 6.0E-34*TEMP3(IL)**(-2.40)*O2*RHO(IL)*RHO(IL)       
!   REACTION 2 O+O3 --> O2+O2
      RATES(IL,2) = 8.0E-12*EXP(-2060*VTEMP(IL))
!   REACTION 3 OSD+M --> O+M N2 & O2 3RD BODY RATES ALLOWED FOR       
      RATES(IL,3) = RHO(IL)*(6.93E-12*EXP(70*VTEMP(IL)) +              &
                               1.677E-11*EXP(110*VTEMP(IL)))
!   REACTION 4 NO+O3 --> NO2+O2
      RATES(IL,4) = 3.0E-12*EXP(-1500*VTEMP(IL))
!   REACTION 5 NO2+O --> NO+O2
      RATES(IL,5) = 5.1E-12*EXP(210*VTEMP(IL))
!   REACTION 6 NO2+O3 --> NO3+O2
      RATES(IL,6) = 1.2E-13*EXP(-2450*VTEMP(IL))
!   REACTION 7 NO2+NO3 + M --> N2O5 + M
      R0 = 2.0E-30*RHO(IL)*TEMP3(IL)**(-4.4)
      R1 = 1.4E-12*TEMP3(IL)**(-0.7)
      RX = ALOG10(R0/R1)
      RATES(IL,7) =(R0/(1.0 + R0/R1))*                                &
                          (0.6**(1.0/(1.0 + (RX*RX))))
!   REACTION 8 N2O5+M --> NO3+NO2+M
      RATES(IL,8) = RATES(IL,7)/(2.7E-27*EXP(11000.0*VTEMP(IL)))
!   REACTION 9 OSD+H2O --> OH+OH
      RATES(IL,9) = 1.63E-10*EXP(60.0*VTEMP(IL))
!   REACTION 10 OH+O3 --> HO2+O2
      RATES(IL,10) = 1.7E-12*EXP(-940.0*VTEMP(IL))
!   REACTION 11 O+HO2 --> OH+O2
      RATES(IL,11) = 3.0E-11*EXP(200.0*VTEMP(IL))
!   REACTION 12 O+OH --> H+O2
      RATES(IL,12) = 2.2E-11*EXP(120.0*VTEMP(IL))
!   REACTION 13 H+O3 --> OH+O2
      RATES(IL,13) = 1.4E-10*EXP(-470.0*VTEMP(IL))
!   REACTION 14 H+O2+M --> HO2+M : MULTIPLIED BY O2 NUMBER DENSITY    
      R0 = 4.4E-32*RHO(IL)*TEMP3(IL)**(-1.3)
      R1 = 4.7E-11*TEMP3(IL)**(-0.2)
      RX = ALOG10(R0/R1)
      RATES(IL,14) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
      RATES(IL,14) = RATES(IL,14)*O2*RHO(IL)
!   REACTION 15 HO2+O3 --> OH+O2+O2
      RATES(IL,15) = 1.0E-14*EXP(-490.0*VTEMP(IL))
!   REACTION 16 OH+HO2 --> H2O+O2
      RATES(IL,16) = 4.8E-11*EXP(250.0*VTEMP(IL))
!   REACTION 17 OH+NO2+M --> HONO2+M
!   INCLUDE ALSO OH+NO2+M --> HOONO + M
      R0 = 1.8E-30*RHO(IL)*TEMP3(IL)**(-3.0)
      R1 = 2.8E-11
      RX = ALOG10(R0/R1)
      RATES(IL,17) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
      R0 = 9.1E-32*RHO(IL)*TEMP3(IL)**(-3.9)
      R1 = 4.2E-11*TEMP3(IL)**(-0.5)
      RX = ALOG10(R0/R1)
      RATES(IL,17) = RATES(IL,17) + (R0/(1.0 + R0/R1))*           &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 18 NO+HO2 --> NO2+OH
      RATES(IL,18) = 3.5E-12*EXP(250.0*VTEMP(IL))
!   REACTION 19 OH+HNO3 --> H2O+NO3 : SPECIAL TREATMENT               
      RX1 = 6.5E-34*RHO(IL)*EXP(1335.0*VTEMP(IL))            
      RX2 = 2.7E-17*EXP(2199.0*VTEMP(IL))            
      RATES(IL,19) =  2.4E-14*EXP(460.0*VTEMP(IL)) +                        &
                      RX1/(1.0 + RX1/RX2)           
!   REACTION 20 CL+O3 --> CLO+O2
      RATES(IL,20) = 2.3E-11*EXP(-200.0*VTEMP(IL))
!   REACTION 21 O+CLO --> CL+O2
      RATES(IL,21) = 2.8E-11*EXP(85.0*VTEMP(IL))
!   REACTION 22 CLO+NO --> NO2+CL
      RATES(IL,22) = 6.4E-12*EXP(290.0*VTEMP(IL))
!   REACTION 23 CLO+NO2+M --> CLONO2+M
      R0 = 1.8E-31*RHO(IL)*TEMP3(IL)**(-3.4)
      R1 = 1.5E-11*TEMP3(IL)**(-1.9)
      RX = ALOG10(R0/R1)
      RATES(IL,23) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 24 O+CLONO2--> CLO+NO3
      RATES(IL,24) = 2.9E-12*EXP(-800.0*VTEMP(IL))
!   REACTION 25 CL+CH4 (+O2) --> HCL+CH3O2
      RATES(IL,25) = 7.3E-12*EXP(-1280.0*VTEMP(IL))
!   REACTION 26 OH+HCL --> H2O+CL
      RATES(IL,26) = 2.6E-12*EXP(-350.0*VTEMP(IL))
!   REACTION 27 CL+HO2 --> HCL+O2
      RATES(IL,27) = 1.8E-11*EXP(170.0*VTEMP(IL))
!   REACTION 28 CLO+HO2 --> HOCL+O2
      RATES(IL,28) = 2.7E-12*EXP(220.0*VTEMP(IL))
!   REACTION 29 CLO+OH --> HO2+CL
      RATES(IL,29) = 7.4E-12*EXP(270.0*VTEMP(IL))
!   REACTION 30 OSD+CH4 --> OH + CH3O2
      RATES(IL,30) = 1.5E-10*EXP(0.0*VTEMP(IL))
!   REACTION 31 OH+CH4 --> H2O + CH3O2
      RATES(IL,31) = 2.45E-12*EXP(-1775.0*VTEMP(IL))
!   REACTION 32 H2CO+OH (+O2) --> H2O + HO2 + CO
      RATES(IL,32) = 5.5E-12*EXP(125.0*VTEMP(IL))
!   REACTION 33 H2CO+CL (+O2) --> HCL+ HO2 + CO        
      RATES(IL,33) = 8.1E-11*EXP(-30.0*VTEMP(IL))
!   REACTION 34 HO2+HO2 --> H2O2+O2
      RATES(IL,34) = 3.5E-13*EXP(430.0*VTEMP(IL)) +                        &
                     1.7E-33*EXP(1000.0*VTEMP(IL))*RHO(IL)                
!   REACTION 35 OH+OH+M --> H2O2+M
      R0 = 6.9E-31*RHO(IL)*TEMP3(IL)**(-1.0)
      R1 = 2.6E-11*TEMP3(IL)**0.0
      RX = ALOG10(R0/R1)
      RATES(IL,35) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 36 OH+H2O2 --> HO2+H2O
      RATES(IL,36) = 2.9E-12*EXP(-160.0*VTEMP(IL))
!   REACTION 37 CLO+CLO+M --> CL2O2 + M
      R0 = 1.6E-32*RHO(IL)*TEMP3(IL)**(-4.5)
      R1 = 2.0E-12*TEMP3(IL)**(-2.4)
      RX = ALOG10(R0/R1)
      RATES(IL,37) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 38 CL2O2+M --> 2CLO + M
      RATES(IL,38) = RATES(IL,37)/(9.3E-28*EXP(8835.0*VTEMP(IL)))  
!   REACTION 39 BR + O3 --> BRO + O2
      RATES(IL,39) = 1.7E-11*EXP(-800.0*VTEMP(IL))
!   REACTION 40 BRO + NO2 + M --> BRONO2 + M
      R0 = 5.2E-31*RHO(IL)*TEMP3(IL)**(-3.2)
      R1 = 6.9E-12*TEMP3(IL)**(-2.9)
      RX = ALOG10(R0/R1)
      RATES(IL,40) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 41 BRO + CLO --> BR + CL + 02
      RATES(IL,41) =  2.3E-12*EXP(260.0*VTEMP(IL))   
!   REACTION 42 BRO + HO2 --> HOBR + O2
      RATES(IL,42) = 4.5E-12*EXP(460.0*VTEMP(IL))
!   REACTION 43 BRO + NO --> BR + NO2
      RATES(IL,43) = 8.8E-12*EXP(260.0*VTEMP(IL))
!   REACTION 44 HOBR + O --> BRO + OH
      RATES(IL,44) = 1.2E-10*EXP(-430.0*VTEMP(IL))
!   REACTION 45 BR + HO2 --> HBR + O2
      RATES(IL,45) = 4.8E-12*EXP(-310.0*VTEMP(IL))
!   REACTION 46 BR + H2CO (+O2) --> HBR + HO2 + CO
      RATES(IL,46) = 1.7E-11*EXP(-800.0*VTEMP(IL))
!   REACTION 47 HBR + OH --> BR + H2O
      RATES(IL,47) = 5.5E-12*EXP(200.0*VTEMP(IL))
!   REACTION 48 BRO + CLO --> BRCL + O2
      RATES(IL,48) = 4.1E-13*EXP(290.0*VTEMP(IL))
!   REACTION 49 HO2 + NO2 + M --> HNO4 + M
      R0 = 2.0E-31*RHO(IL)*TEMP3(IL)**(-3.4)
      R1 = 2.9E-12*TEMP3(IL)**(-1.1)
      RX = ALOG10(R0/R1)
      RATES(IL,49) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 50 HNO4 + OH --> H2O + NO2 + O2
      RATES(IL,50) = 1.3E-12*EXP(380.0*VTEMP(IL))
!
! REACTIONS ADDED FOR TROPOSPHERIC CHEMISTRY
!
!   REACTION 51 OH + OH --> H2O + O
      RATES(IL,51) = 1.8E-12
!   REACTION 52 CO + OH (+O2) --> HO2 +CO2
      R0 = 5.9E-33*RHO(IL)*TEMP3(IL)**(-1.4)
      R1 = 1.1E-12*TEMP3(IL)**(1.3)
      RX = ALOG10(R0/R1)
      RATES(IL,52) =    (R0/(1.0 + R0/R1))*                           &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
      R0 = 1.5E-13*TEMP3(IL)**0.6
      R1 = 2.1E9*TEMP3(IL)**(6.1)/RHO(IL)       
      RX = ALOG10(R0/R1)       
      RATES(IL,52) = RATES(IL,52)+ (R0/(1.0 + R0/R1))*                 &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 53 NO3 + NO --> 2NO2
      RATES(IL,53) = 1.5E-11*EXP(170.0*VTEMP(IL))
!   REACTION 54 H2 + OH (+O2) --> HO2 + H2O
      RATES(IL,54) = 5.5E-12*EXP(-2000.0*VTEMP(IL))
!   REACTION 55 CH3O2 + NO (+O2) --> H2CO + NO2 + HO2
      RATES(IL,55) = 2.8E-12*EXP(300.0*VTEMP(IL))
!   REACTION 56 CH3O2 + HO2 --> CH3OOH + O2
      RATES(IL,56) = 4.1E-13*EXP(750.0*VTEMP(IL))
!   REACTION 57 CH3OOH + OH --> CH3O2 + H2O 
      RATES(IL,57) = 3.8E-12*EXP(200.0*VTEMP(IL))
!   REACTION 58 NO2 + NO3 --> NO + NO2 + O2
      RATES(IL,58) = 4.5E-14*EXP(-1260.0*VTEMP(IL))
!   REACTION 59 NO3 + NO3 --> 2NO2 + O2 
      RATES(IL,59) = 8.5E-13*EXP(-2450.0*VTEMP(IL))
!   REACTION 60 OH + HONO --> H2O + NO2
      RATES(IL,60) = 1.8E-11*EXP(-390.0*VTEMP(IL))
!   REACTION 61 CH3O2 + CH3O2 --> 2H2CO + 2OH
      RATES(IL,61) = 9.5E-14*EXP(390.0*VTEMP(IL))
!   REACTION 62 NO3 + HCHO (+O2) --> HNO3 + CO + HO2
      RATES(IL,62) = 5.8E-16
!   REACTION 63 HO2NO2 + M --> HO2 + NO2
      RATES(IL,63) = RATES(IL,49)*4.762E26*EXP(-10900*VTEMP(IL))
!   REACTION 64 NO + O + M --> NO2 + M
      R0 = 9.0E-31*RHO(IL)*TEMP3(IL)**(-1.5)
      R1 = 3.0E-11
      RX = ALOG10(R0/R1)
      RATES(IL,64) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 65 OH + NO + M --> HONO + M
      R0 = 7.0E-31*RHO(IL)*TEMP3(IL)**(-2.6)
      R1 = 3.6E-11*TEMP3(IL)**(-0.1)
      RX = ALOG10(R0/R1)
      RATES(IL,65) =(R0/(1.0 + R0/R1))*                               &
                          (0.6**(1.0/(1.0 + (RX*RX)))) 
!   REACTION 66 N2O + OSD --> 2NO
      RATES(IL,66) = 6.7E-11*EXP(20.0*VTEMP(IL))
!   REACTION 67 N + NO --> N2 + O
      RATES(IL,67) = 2.1E-11*EXP(100.0*VTEMP(IL))
!   REACTION 68 N + O2 --> NO + O
      RATES(IL,68) = 1.5E-11*EXP(-3600.0*VTEMP(IL))
!   REACTION 69 CLO+OH --> --> HCL+O2
      RATES(IL,69) = 6.0E-13*EXP(-230.0*VTEMP(IL))
!   REACTION 70 N2O + OSD --> N2 + O2
      RATES(IL,70) = 4.7E-11*EXP(20.0*VTEMP(IL))
!   REACTION 71 N + NO2 --> N2O + O
      RATES(IL,71) = 5.8E-12*EXP(220.0*VTEMP(IL))
  ENDDO
  DO IC=IHET,IHET+NHET-1
    DO IL=1,MERIDS
      RATES(IL,IC) = 0.0                                            
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE REACTION_RATES

SUBROUTINE PSC(TEMP,PRESS,RHO,DYY,H2O,H2SO4,ANAT,AICE,COND,TICE,  &
           WH2SO4,AM,AW,ALIQ,RMEAN,ASAT,RNAT,RICE,MERIDS,NCHEM,MYPE)
!---------------------------------------------------------------------------
!  CALCULATES IMPORTANT TERMS FOR BINARY AEROSOL/SAT/NAT/ICE SCHEME.
! THE MAIN OUTPUT IS THE SURFACE AREA OR PARTICLE RADIUS
!----------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, intent(in)                           :: MERIDS,NCHEM,MYPE
REAL, intent(in)                              :: ANAT, AICE
REAL, intent(in),    dimension(MERIDS)        :: PRESS, RHO, TEMP, H2SO4

REAL, intent(inout), dimension(MERIDS)        :: H2O, TICE, WH2SO4, AM, AW, &
                                                 ALIQ, RMEAN, ASAT, RNAT, RICE
REAL, intent(inout), dimension(MERIDS,3)      :: COND
REAL, intent(inout), dimension(MERIDS, NCHEM) :: DYY

INTEGER IL
REAL DENS(MERIDS),VHET(MERIDS,4),WF(MERIDS,2),PH2O(MERIDS)
REAL AVGDR,RR,PI,ADROP,SIGMAL,TSAT,ABT,AMT,BMT,PX,C2,   &
     A1,A2,A3,A4,C3,C4,CLIMIT,P0H2O,Y1,Y2,T1,RMODE,RMODESAT,SANAT
!----------------------------------------------------------------------------
!
! AVGDR IS THE RECIPROCAL OF THE AVOGADRO CONSTANT
! RR IS THE GAS CONSTANT
!
!----------------------------------------------------------------------------
  DATA AVGDR,RR/1.66056547E-24, 8.3144/,PI/3.1415926536/
!----------------------------------------------------------------------------
!
!   ADROP, ANAT AND AICE ARE THE (FIXED) NUMBER OF DROPS OR PARTICLES PER CC
!   IN THE POLAR STRATOSPHERC CLOUDS. A LOG NORMAL SIZE DISTRIBUTION IS ASSUMED
!   WITH SIGMA = 1.8 GIVING LOG(SIGMA)**2 = 0.3455  
!
!----------------------------------------------------------------------------
  DATA ADROP/10.0/,SIGMAL/0.34549316/      
!          
    DO 40 IL=1,MERIDS
!----------------------------------------------------------------------------
!
! COMPUTE SOLID CONDENSED MATTER: COND(IL,1:3) -- SAT, NAT, ICE
!               
!----------------------------------------------------------------------------
      PH2O(IL) = H2O(IL)*PRESS(IL)/(RHO(IL)*101325.0)
      COND(IL,1) = 0.0
      COND(IL,2) = 0.0
      COND(IL,3) = 0.0
      TSAT = 3236.0/(11.502 - ALOG10(PH2O(IL)*760.0))           
      IF(TEMP(IL).LT.TSAT) COND(IL,1) = H2SO4(IL)        
      ABT = 39.1104 - 11397.0/TEMP(IL) + 9.179E-3*TEMP(IL)                
      AMT = -2.7836 - 8.8E-4*TEMP(IL)
      BMT = -2.1249 + ALOG10(PH2O(IL)*101325.0)
      PX = AMT*BMT + ABT
      C2 = DYY(IL,1) - (RHO(IL)*100.0/PRESS(IL))*10.0**(PX)
      IF(C2.GT.0.0) COND(IL,2) = C2
      TICE(IL) = 2668.70/(10.4310 - ALOG10(760.0*PH2O(IL)))               
      A1 = 7.5502 - 2668.7/TEMP(IL)
      C3 = H2O(IL) - (RHO(IL)*101325.0/PRESS(IL))*10.0**(A1)          
      C4 = (RHO(IL)*101325.0/PRESS(IL))*10.0**(A1)          
      IF(C3.GT.0.0) COND(IL,3) = C3
      DYY(IL,1) = DYY(IL,1) - COND(IL,2)
      DYY(IL,15) = DYY(IL,15) - COND(IL,2) 
      CLIMIT = RHO(IL)*1.E-15
      IF(DYY(IL,15).LT.CLIMIT) DYY(IL,15) = CLIMIT
      H2O(IL) = H2O(IL) - COND(IL,3) 
 40 ENDDO!CONTINUE
!-------------------------------------------------------------------------
!
!  COMPUTE WEIGHT % H2SO4. FROM TABAZADEHET AL., GRL, 24, 1931-1934, 1997.
!  TABLE A1 OF SHIAET AL. JGR, 106, 24,529-24,274, 2001
!
!-------------------------------------------------------------------------
    DO 50 IL = 1,MERIDS
      P0H2O =EXP(18.452406985 - 3505.1578807/TEMP(IL) -                &
             330918.55082/(TEMP(IL)**2) + 12725068.262/(TEMP(IL)**3))    
      AW(IL) = H2O(IL)*PRESS(IL)*0.01/(RHO(IL)*P0H2O)
      IF(AW(IL).LE.0.05) THEN
        Y1 = 12.37208932*AW(IL)**(-0.16125516114) - 30.490657554*AW(IL)  &
            -  2.1133114241
        Y2 = 13.455394705*AW(IL)**(-0.1921312255) - 34.285174607*AW(IL)        &
            - 1.7620073078
      ELSEIF(AW(IL).LT.0.85.AND.AW(IL).GT.0.05) THEN
        Y1 = 11.820654354*AW(IL)**(-0.20786404244) - 4.807306373*AW(IL)   &
            -  5.1727540348
        Y2 = 12.891938068*AW(IL)**(-0.23233847708) - 6.4261237757*AW(IL)   &
            - 4.9005471319
      ELSE
        Y1 = -180.06541028*AW(IL)**(-0.38601102592)                     &
            - 93.317846778*AW(IL) + 273.88132245
        Y2 = -176.95814097*AW(IL)**(-0.36257048154)                     &
             - 90.469744201*AW(IL) + 267.45509988
      ENDIF              
      AM(IL) = Y1 + (TEMP(IL) - 190.0)*(Y2 - Y1)/70.0
      WH2SO4(IL) = 9800.0*AM(IL)/(98.0*AM(IL) + 1000.0)
      WF(IL,1) = 0.01*WH2SO4(IL) 
      WF(IL,2) = 0.0
 50 ENDDO!CONTINUE
!---------------------------------------------------------------------------
!
!  COMPUTE DENSITY OF BINARY AEROSOL        
!
!---------------------------------------------------------------------------
    CALL DENSITY(WF,TEMP,DENS,MERIDS)
!---------------------------------------------------------------------------
!
!  COMPUTE VOLUME OF BINARY AEROSOL/SAT/NAT/ICE 
!  1.6, 1.35 and 0.928  are the densities of SAT, NAT and ICE     
!
!---------------------------------------------------------------------------
    DO 100 IL=1,MERIDS
      T1 = H2SO4(IL)*PRESS(IL)/(RHO(IL)*TEMP(IL)*RR)
      VHET(IL,1) = T1*98.076E-6/(WF(IL,1)*DENS(IL))             
      VHET(IL,2) = COND(IL,1)*170.1*AVGDR/1.6
      VHET(IL,3) = COND(IL,2)*117.1*AVGDR/1.35
      VHET(IL,4) = COND(IL,3)*18.02*AVGDR/0.928
100 ENDDO
!---------------------------------------------------------------------------
!
!  COMPUTE PARTICLE PARAMETERS FROM WHICH THE HETEROGENEOUS REACTION RATES
!  ARE DETERMINED; ASSUME SURFACE AREA FROM SAT IS LIMITED BY NAT AMOUNT   
!
!---------------------------------------------------------------------------
    A1 =EXP(-4.5*SIGMAL)
    A2 =EXP(0.5*SIGMAL)
    A3 =EXP(2.0*SIGMAL)
    A4 = 1.33333333*PI*ADROP
    DO 150 IL=1,MERIDS
      RMODE = (VHET(IL,1)*A1/A4)**0.33333333
      RMEAN(IL) = RMODE*A2
      ALIQ(IL) = 3.0*A4*(RMODE**2)*A3
      IF(RMEAN(IL).LT.1.0E-12) RMEAN(IL) = 1.0E-12
      RMODESAT = (VHET(IL,2)*A1/A4)**0.33333333
      ASAT(IL) = 3.0*A4*(RMODESAT**2)*A3
      RNAT(IL) = (VHET(IL,3)/(1.33333333*PI*ANAT))**0.33333333
      RICE(IL) = (VHET(IL,4)/(1.33333333*PI*AICE))**0.33333333
      SANAT = 4.0*PI*ANAT*RNAT(IL)**2
      ASAT(IL) = ASAT(IL) - SANAT
      ALIQ(IL) = ALIQ(IL) - SANAT
      IF(ASAT(IL).LT.0.0) ASAT(IL) = 0.0
      IF(ALIQ(IL).LT.0.0) ALIQ(IL) = 0.0
150 ENDDO
    RETURN
END SUBROUTINE PSC
!
SUBROUTINE GAMMA_AER(TEMP,PRESS,RHO,DYY,H2O,WH2SO4,AM,AW,RMEAN,   &
                     GAMMA,RTT,MERIDS,NCHEM,NHET,MYPE)
!-------------------------------------------------------------------------
!
! SUBROUTINE TO CALCULATE REACTION PROBABILITIES ON SULPHATE AEROSOL
! BASED ON JPL'03 RECOMMENDATION
!  
!-------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, intent(in) ::  MERIDS,NCHEM,NHET,MYPE
  REAL, intent(in) :: TEMP(MERIDS),PRESS(MERIDS),RHO(MERIDS),DYY(MERIDS,NCHEM), &
        H2O(MERIDS),WH2SO4(MERIDS),AW(MERIDS),AM(MERIDS),           &
        RMEAN(MERIDS)
  REAL, intent(inout) :: RTT(MERIDS), GAMMA(MERIDS,NHET)
!
  INTEGER  IL
  REAL  AMH2SO4(MERIDS),XH2SO4(MERIDS),VISC(MERIDS),AACID(MERIDS),  &
        TEMP2(MERIDS)
!-------------------------------------------------------------------------
!  
! CALCULATE H2SO4 MOLARITY (AMH2SO4), MOLE FRACTION (XH2SO4), 
!  VISCOSITY (VISC) AND ACID ACTIVITY (AACID)
!  TABLE A2 OF SHIAET AL. JGR, 106, 24,529-24,274, 2001.
! 
!-------------------------------------------------------------------------
  REAL  T2,Z1,Z2,Z3,RHOX,AA,X,T1,T3,AKH,AKH2O,AKHYDR,DIFF,SCLONO2,  &
         CCLONO2,GAMMAB1,HHCL,AKHCL,Q1,RQ,A1,FCLONO2,GAMMARXN,      &
         GAMMABHCL,GAMMAS,FHCL,GAMMASP,GAMMABHCLP,GAMMAB,GCLONO2,   &
         SHOCL,HHOCL,FHOCL,WT,AK0,AK1,AK2,T0,HCLONO2,   &
         AMHCL,AKHOCL,CHOCL
!
!  The parameterisations used here break down below about 185K, so the
!  temperature is here limited to 185K and above (TEMP2).
!
    DO IL = 1,MERIDS
      TEMP2(IL) = TEMP(IL)
      IF(TEMP2(IL).LT.185.0) TEMP2(IL) = 185.0
      RTT(IL) = SQRT(TEMP(IL))
    ENDDO
    DO 20 IL = 1,MERIDS
      T2 = TEMP2(IL)**2
      Z1 = 0.12364 - 5.6E-7*T2
      Z2 = -0.02954 + 1.814E-7*T2
      Z3 = 2.343E-3 - 1.487E-6*TEMP2(IL) - 1.324E-8*T2
      RHOX = 1.0 + Z1*AM(IL) + Z2*AM(IL)**1.5 + Z3*AM(IL)**2
      AMH2SO4(IL) = RHOX*WH2SO4(IL)/9.8
      XH2SO4(IL) = WH2SO4(IL)/                                          &
         (WH2SO4(IL) + (100.0 - WH2SO4(IL))*98.0/18.0)
      AA = 169.5 + 5.18*WH2SO4(IL) - 0.0825*WH2SO4(IL)**2 +             &
           3.27E-3*WH2SO4(IL)**3 
      T0 = 144.11 + 0.166*WH2SO4(IL) - 0.015*WH2SO4(IL)**2 +            &
           2.18E-4*WH2SO4(IL)**3
      X = TEMP2(IL)**(-1.43) 
      VISC(IL) = AA*X*EXP(448.0/(TEMP2(IL) - T0))
      T1 = 60.51 - 0.095*WH2SO4(IL) + 0.0077*WH2SO4(IL)**2              &
           - 1.61E-5*WH2SO4(IL)**3
      T2 =  (-805.89 + 253.05*WH2SO4(IL)**0.076)/RTT(IL)  
      T3 =   (1.76 + 2.52E-4*WH2SO4(IL)**2)*RTT(IL)
      AACID(IL) =EXP(T1 + T2 - T3)
 20 ENDDO!CONTINUE
!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR CLONO2 + H2O AND CLONO2 + HCL AND 
!  HENRY'S LAW COEFFICIENTS.
!  TABLE A3 OF SHIAET AL. JGR, 106, 24,529-24,274, 2001.
!
!  The following formulation for the water activity is from Shi et al.,
!  but the differences between their parameterisation and that calculated
!  using the actua model H2O is not large.
!
!      awx = exp((-69.775*xh2so4(il) - 18253.7*xh2so4(il)**2 + 
!     2     31072.2*xh2so4(il)**3 - 25668.8*xh2so4(il)**4)*
!     3     (1.0/temp(il) - 26.9033/(temp(il)**2)))      
!      AKHYDR = AWx*(AKH2O + AKH*AACID(IL))
!
!-------------------------------------------------------------------------
    DO 30 IL = 1,MERIDS
      AKH = 1.22E12*EXP(-6200.0/TEMP2(IL))
      AKH2O = 1.95E10*EXP(-2800.0/TEMP2(IL))     
      AKHYDR = AW(IL)*(AKH2O + AKH*AACID(IL))
      DIFF = 5.0E-8*TEMP2(IL)/VISC(IL)
      SCLONO2 = 0.306 + 24.0/TEMP2(IL)
      HCLONO2 = 1.6E-6*EXP(4710.0/TEMP2(IL) - SCLONO2*AMH2SO4(IL))
      CCLONO2 = 1474.*TEMP2(IL)**0.5
      GAMMAB1 = (4.0*HCLONO2*0.082*TEMP2(IL)/CCLONO2)*                  &
                (DIFF*AKHYDR)**0.5
!
      HHCL = (0.094 - 0.61*XH2SO4(IL) + 1.2*XH2SO4(IL)**2)*             &
        EXP(-8.68 + (8515. - 10718.*XH2SO4(IL)**0.7)/TEMP2(IL))
      AMHCL = HHCL*DYY(IL,4)*PRESS(IL)/(RHO(IL)*101325.0)
      AKHCL = 7.9E11*AACID(IL)*DIFF*AMHCL
      Q1 = (DIFF/(AKHYDR + AKHCL))**0.5
      RQ = RMEAN(IL)/Q1
      A1 = RQ + 0.312*RQ**2
      FCLONO2 = A1/(3.0 + A1)
      GAMMARXN = FCLONO2*GAMMAB1*(1.0 + AKHCL/AKHYDR)**0.5
      GAMMABHCL = GAMMARXN*AKHCL/(AKHCL + AKHYDR)
      GAMMAS = 66.12*EXP(-1374./TEMP2(IL))*HCLONO2*AMHCL
      IF(DYY(IL,4).NE.0.0) THEN
        FHCL = 1.0/(1.0 + 0.612*(GAMMAS + GAMMABHCL)*DYY(IL,6)/DYY(IL,4))
      ELSE
        FHCL = 0.0
      ENDIF
      GAMMASP = FHCL*GAMMAS
      GAMMABHCLP = FHCL*GAMMABHCL     
      GAMMAB = GAMMABHCLP + GAMMARXN*AKHYDR/(AKHCL + AKHYDR)
      GCLONO2 = 1.0/(1.0 + 1.0/(GAMMASP + GAMMAB))
      GAMMA(IL,5) = GCLONO2*(GAMMASP + GAMMABHCLP)/                     &
            (GAMMASP + GAMMAB)
      GAMMA(IL,4) = GCLONO2 - GAMMA(IL,5)
!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR HOCL + HCL AND HENRY'S LAW COEFFICIENTS.
!  TABLE A4 OF SHIAET AL. JGR, 106, 24,529-24,274, 2001.
!
!-------------------------------------------------------------------------
      SHOCL = 0.0776 + 59.18/TEMP2(IL)
      HHOCL = 1.91E-6*EXP(5862.4/TEMP2(IL) - SHOCL*AMH2SO4(IL))
      DIFF = 6.4E-8 *TEMP2(IL)/VISC(IL) 
      AKHOCL = 1.25E9*AACID(IL)*DIFF*AMHCL
      CHOCL = 2009.*RTT(IL)
      Q1 = (DIFF/AKHOCL)**0.5
      RQ = RMEAN(IL)/Q1
      A1 = RQ + 0.312*RQ**2
      FHOCL = A1/(3.0 + A1)
      GAMMARXN = (FHOCL*4.0*HHOCL*0.082*TEMP2(IL)/CHOCL)*               &
                 (DIFF*AKHOCL)**0.5
      GAMMA(IL,1) = 1.0/(1.0 + 1.0/(GAMMARXN*FHCL))
 30 ENDDO!CONTINUE
!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR N2O5 + H2O 
!  ROBINSONET AL. JGR, 102, 3583-3601, 1997.
!
!-------------------------------------------------------------------------
    DO 40 IL = 1,MERIDS
      WT = WH2SO4(IL)
      IF(WH2SO4(IL).GT.80.0) WT = 80.0
      AK0 = -25.5265 - 0.133188*WT + 0.00930846*WT**2 -                 &
              9.0194E-5*WT**3  
      AK1 = 9283.76 + 115.345*WT - 5.19258*WT**2 +                      &
           0.0483464*WT**3  
      AK2 = -851801. - 22191.2*WT + 766.916*WT**2 -                     &
           6.85427*WT**3
      GAMMA(IL,3) =EXP(AK0 + AK1/TEMP2(IL) + AK2/(TEMP2(IL)**2))  
 40 ENDDO!CONTINUE
!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR 
!       N2O5 + HCL
!       HOBR + HCL
!       HOCL + HBR 
!  NO RECOMMENDATION IN JPL '03, ASSUMED ZERO
!
!-------------------------------------------------------------------------
    DO IL = 1,MERIDS
      GAMMA(IL,2) = 0.0
      GAMMA(IL,6) = 0.0
      GAMMA(IL,7) = 0.0
    ENDDO
!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR HOBR + HBR 
!  ABBATT, JGR, 100, 14009-14017, 1995. 
!
!-------------------------------------------------------------------------
    DO IL = 1,MERIDS
      GAMMA(IL,8) = 0.25
    ENDDO
!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR BRONO2 + H2O 
!  USE JPL '03 RECOMMENDATION (HANSON PERS. COMM.)
!
!-------------------------------------------------------------------------
    DO IL = 1,MERIDS
      GAMMA(IL,9) = 1.0/(1.2422 + 1.0/(0.114 +EXP(29.24 -              &
                    0.396*WH2SO4(IL))))
    ENDDO
    RETURN
END SUBROUTINE GAMMA_AER
!
SUBROUTINE HETRATES(TEMP,PRESS,RHO,DYY,H2O,TICE,ANAT,AICE,ALIQ,        &
        RMEAN,RNAT,RICE,GAMMA,RTT,RATES,IHET,NHET,MERIDS,NCHEM,MYPE)

IMPLICIT NONE
!------------------------------------------------------------------------
!
!  This subroutine computes the equivalent 2nd order reaction rates for 
!  the heterogeneous reactions on aerosol, nat and ice.
!
!------------------------------------------------------------------------
      
  INTEGER, intent(in) :: IHET, NHET, MERIDS, NCHEM, MYPE
  REAL, intent(in) ::  TEMP(MERIDS),PRESS(MERIDS),RHO(MERIDS), &
       DYY(MERIDS,NCHEM), H2O(MERIDS),TICE(MERIDS),            &
       ALIQ(MERIDS),RMEAN(MERIDS),RNAT(MERIDS),RICE(MERIDS),   &
       GAMMA(MERIDS,NHET), RTT(MERIDS)
  REAL, intent(inout) :: RATES(MERIDS,IHET+NHET-1)
!
  INTEGER INN(9),IC,IL
  REAL AMW(9),GNAT(9),GICE(9),                          &
       CHEMC(MERIDS,NHET),CONST(NHET),                              &
       G2NAT(MERIDS,NHET),DELT(MERIDS),SICE(MERIDS)
!
  REAL  ANAT,AICE,PI,ANUM,ADEN
  DATA AMW/52.45, 108.00, 108.00, 97.45, 97.45, 96.91,              &
           52.45, 96.91, 141.91/               
!------------------------------------------------------------------------
!  AMW = MOLECULAR WEIGHT OF GAS PHASE SPECIES
!------------------------------------------------------------------------
  DATA PI/3.141592653589793/    
  DATA GNAT/0.1, 3.0E-3, 4.0E-4, 4.0E-3, 0.2,                    &
            0.0, 0.0, 0.0, 0.0/
  DATA GICE/0.2, 0.03, 0.02, 0.3, 0.3, 0.3,                         &
            0.03, 0.1, 0.3/
  DATA INN/4,4,0,0,4,4,11,11,0/
!------------------------------------------------------------------------
!  INN = INDEX NUMBER OF LIQUID/SOLID PHASE SPECIES (0= H2O)
!-------------------------------------------------------------------------
!     REACTION 70 HOCL + HCL --> H2O + CL2 (HETEROGENEOUS)              
!     REACTION 71 N2O5 + HCL --> HNO3 + CLNO2 (HETEROGENEOUS)           
!     REACTION 72 N2O5+H2O --> 2HNO3 (HETEROGENEOUS)
!     REACTION 73 CLONO2+H2O --> HOCL+HNO3 (HETEROGENEOUS)              
!     REACTION 74 CLONO2+HCL --> CL2+HNO3 (HETEROGENEOUS)               
!     REACTION 75 HOBR + HCL --> BRCL + H2O (HETEROGENEOUS)             
!     REACTION 76 HOCL + HBR --> BRCL + H2O (HETEROGENEOUS)             
!     REACTION 77 HOBR + HBR --> 2BR + H2O (HETEROGENEOUS)              
!     REACTION 78 BRONO2 + H2O --> HOBR + HNO3 (HETEROGENEOUS)
!     (THE FIRST MOLECULE ON THE LEFT HAND SIDE IS IN THE GAS PHASE,
!     THE SECOND MOLECULE IS IN THE LIQUID/SOLID PHASE)         
!
!-------------------------------------------------------------------------
!
!    aliq is the liquid surface area. const is sqrt(8R/(pi*mw)) for the
!    gaseous phase species, with the mean molecular
!    speed equal to const*sqrt(temp)
!
!-------------------------------------------------------------------------
    DO 60 IC = 1,NHET
      CONST(IC) = SQRT(8.0*8.3144E7/(PI*AMW(IC)))
      DO 50 IL = 1,MERIDS
        IF(INN(IC).EQ.0) THEN
          CHEMC(IL,IC) = H2O(IL)
        ELSE
          CHEMC(IL,IC) = DYY(IL,INN(IC))
        ENDIF
        RATES(IL,IHET-1+IC) =  0.0 
 50   ENDDO
 60 ENDDO
!-------------------------------------------------------------------------
!
!    Reactions IHET to IHET+NHET-1 on NAT
!
!-------------------------------------------------------------------------
    DO IC = 1,NHET
      DO IL = 1,MERIDS
        G2NAT(IL,IC) = GNAT(IC)
      ENDDO
    ENDDO
    DO 110 IL=1,MERIDS
      DELT(IL) = TEMP(IL) - TICE(IL)
      SICE(IL) =  10**(2668.70*(1.0/TEMP(IL) - 1.0/TICE(IL)))      
      IF (SICE(IL) .GT. 3.0) SICE(IL) = 3.0
      G2NAT(IL,4) =EXP(-9.03 + 2.81*SICE(IL))
      G2NAT(IL,5) = 1.0/(4.3478 + 1.4241*EXP(0.518*DELT(IL)))             
110 ENDDO
!
    DO IC=1,NHET
      DO IL=1,MERIDS
        ANUM = PI*CONST(IC)*RTT(IL)*G2NAT(IL,IC)*ANAT*RNAT(IL)**2
        ADEN = CHEMC(IL,IC)
        IF(ADEN.GT.0.0.AND.ANUM.GT.0.0)                               &
            RATES(IL,IHET-1+IC) =  ANUM/ADEN
      ENDDO
    ENDDO
!-------------------------------------------------------------------------
!
!    Reactions IHET to IHET+NHET-1 on ICE
!
!------------------------------------------------------------------------ 
    DO IC=1,NHET
      DO IL=1,MERIDS
        ANUM = PI*CONST(IC)*RTT(IL)*GICE(IC)*AICE*RICE(IL)**2  
        ADEN = CHEMC(IL,IC)                
        IF(ADEN.GT.0.0.AND.ANUM.GT.0.0)                               &
           RATES(IL,IHET-1+IC) = ANUM/ADEN 
      ENDDO
    ENDDO
!-------------------------------------------------------------------------
!
!    Reactions IHET to IHET+NHET-1 on AEROSOL
!    aliq is the liquid surface area. const is sqrt(8R/(pi*mw)) for the
!    gaseous phase species, with the mean molecular speed equal to 
!        const*sqrt(temp)
!
!-------------------------------------------------------------------------
    DO IC = 1,NHET
      DO IL = 1,MERIDS
        IF(CHEMC(IL,IC).GT.0.0)                                       &
           RATES(IL,IHET+IC-1) =  RATES(IL,IHET+IC-1) +               &
           0.25*GAMMA(IL,IC)*CONST(IC)*RTT(IL)*ALIQ(IL)/CHEMC(IL,IC)
      ENDDO
    ENDDO
!    
    RETURN
END SUBROUTINE HETRATES
!
SUBROUTINE DENSITY(WF,T,DENS,MERIDS)
IMPLICIT NONE
!
!    Density of ternary solution in g cm-3
!
  INTEGER, intent(in)    :: MERIDS
  REAL,    intent(in )   :: WF(MERIDS,2), T(MERIDS)
  REAL,    intent(inout) :: DENS(MERIDS)

  INTEGER IL
  REAL X(22),AMR(3)
  REAL W,WH,T2,V1,A1,A2,VS,VN,VMCAL
  DATA X/2.393284E-02,-4.359335E-05,7.961181E-08,0.0,-0.198716351,   &
         1.39564574E-03,-2.020633E-06,0.51684706,-3.0539E-03,        &
         4.505475E-06,-0.30119511,1.840408E-03,-2.7221253742E-06,    &
        -0.11331674116,8.47763E-04,-1.22336185E-06,0.3455282,        &
        -2.2111E-03,3.503768245E-06,-0.2315332,1.60074E-03,          &
        -2.5827835E-06/
  DATA AMR/0.05550622,0.01019576,0.01586899/
    DO 100 IL=1,MERIDS
      W = WF(IL,1) + WF(IL,2)
      WH = 1.0 - W
      T2 = T(IL)**2
      V1 = X(1) + X(2)*T(IL) + X(3)*T2 + X(4)*T2*T(IL)
      A1 = X(8) + X(9)*T(IL) + X(10)*T2
      A2 = X(11) + X(12)*T(IL) + X(13)*T2
      VS = X(5) + X(6)*T(IL) + X(7)*T2 + A1*W + A2*W**2
      A1 = X(17) + X(18)*T(IL) + X(19)*T2
      A2 = X(20) + X(21)*T(IL) + X(22)*T2
      VN = X(14) + X(15)*T(IL) + X(16)*T2 + A1*W + A2*W**2
      VMCAL = WH*V1*AMR(1) + VS*WF(IL,1)*AMR(2) + VN*WF(IL,2)*AMR(3)
      DENS(IL) = 1.0E-3/VMCAL
100 ENDDO
!
  RETURN
END SUBROUTINE DENSITY
     
SUBROUTINE PHOTO_RATES (VTEMP,OZONX,COSPHI,KL,PHOTO,OZON,COSP,    &
                        COSPHC,MERIDS,IRX,IRX4,FOTO,MYPE)
!---------------------------------------------------------------------- 
!   THIS SUBROUTINE CALCULATES THE PHOTOLYSIS RATES USING A LOOK-UP     
!   TABLE. THE RATES ARE PUT INTO FOTO AND ARE AS FOLLOWS               
!   1.  O2 + HV --> 2O
!   2.  O3 + HV --> OSD + O2
!   3.  O3 + HV --> O3P + O2
!   4.  NO2 + HV --> NO + O
!   5.  HNO3 + HV --> OH + NO2 (300K) (LATER OVERWRITTEN WITH VALUE     
!   6.  NO3 + HV --> NO2 + O
!   7.  NO3 + HV --> NO + O2
!   8.  N2O5 + HV --> NO2 + NO3  (300K)  (LATER OVERWRITTEN WITH VALUE  
!   9.  N2O5 + HV --> NO2 + NO3  (250K)               AT TEMP)          
!  10.  N2O5 + HV --> NO2 + NO3  (200K)
!  11.  H2O2 + HV --> 2OH
!  12.  CLONO2 + HV --> CL + NO3 (300K) (LATER OVERWRITTEN WITH VALUE   
!  13.  H2CO + HV --> H + HCO --> 2H + CO
!  14.  H2CO + HV --> H2 + CO
!  15.  HOCL + HV --> OH + CL
!  16.  CL2O2 + HV --> CLO2 + CL
!  17.  BRONO2 + HV --> BR + NO3
!  18.  HOBR + HV --> OH + BR
!  19.  BRO + HV --> BR + O
!  20.  HNO3 + HV --> OH + NO2 (250K)
!  21.  HNO3 + HV --> OH + NO2 (200K)
!  22.  CLONO2 + HV --> CL + NO3 (250K)
!  23.  CLONO2 + HV --> CL + NO3 (200K)
!  24.  BRCL + HV --> BR + CL
!  25.  HNO4 + HV --> HO2 + NO2
!  26.  HONO + HV --> OH + NO
!  27.  CH3OOH + HV --> H2CO + HO2 + OH
!  28.  N2O + HV --> N2 + OSD  (300K)  (LATER OVERWRITTEN WITH VALUE  
!  29.  N2O + HV --> N2 + OSD  (250K)               AT TEMP)          
!  30.  N2O + HV --> N2 + OSD  (200K)
!  31.  NO + HV --> N + O
!  32.  H2O + HV --> H + OH
!  33.  HCL + HV --> H + CL
!---------------------------------------------------------------------- 
IMPLICIT NONE
  INTEGER, intent(in) :: MERIDS,IRX,IRX4,MYPE, KL
  REAL, intent(in) :: PHOTO(IRX4,14,11,48), OZON(11,48),      &
       VTEMP(MERIDS),OZONX(MERIDS),COSPHI(MERIDS)
  REAL, intent(inout) ::   COSP(14), FOTO(MERIDS,IRX4), COSPHC(48)


  REAL  X1(MERIDS),X2(MERIDS),Y1(MERIDS),Y2(MERIDS),Z1(MERIDS),       &
        Z2(MERIDS),OZON1(MERIDS),OZON2(MERIDS),           &
        COSP1(MERIDS),COSP2(MERIDS)  
  INTEGER IL,ICOX,IOX,IR,IM
  INTEGER IOZ(MERIDS),ICOZ(MERIDS)    
  LOGICAL LCOZ(MERIDS)
!
!      REAL  ALG,CX1,AX1,BX1,CX2,AX1,BX2,AX3,BX3,XX,ZZ,P1,P2,XX2,AX2,CX3  
  REAL  ALG,CX1,AX1,BX1,CX2,AX2,BX2,CX3,AX3,BX3,XX,ZZ,P1,P2,XX2  
!    
!
     IF(KL.EQ.48) COSPHC(KL) = COSPHC(2)
     COSP(1) = -COSPHC(KL)*0.5
     COSP(2) = -COSPHC(KL)*0.25
     DO 110 IL = 1,MERIDS
       LCOZ(IL) = COSPHI(IL).GT.(-COSPHC(KL))
 110 ENDDO
     DO 120 IL = 1,MERIDS
       ICOZ(IL) = 1
       COSP1(IL) = COSP(1)
       COSP2(IL) = COSP(2)
 120 ENDDO
     DO 140 ICOX =  2,13
       DO 130 IL = 1,MERIDS
         IF(COSPHI(IL).GT.COSP(ICOX))  THEN
           ICOZ(IL) = ICOX
           COSP1(IL) = COSP(ICOX)
           COSP2(IL) = COSP(ICOX+1)
         ENDIF
 130   ENDDO
 140 ENDDO
!
!  INTERPOLATE LINEARLY IN COSPHI AND LINEARLY IN
!  OZONE AMOUNT FOR THE LOG OF THE PHOTODISSOCIATION RATE.              
!
     DO 150 IL = 1,MERIDS
       IOZ(IL) = 1
       OZON1(IL) = OZON(1,KL)
       OZON2(IL) = OZON(2,KL)
 150 ENDDO
     DO 170 IOX = 2,10
       DO 160 IL = 1,MERIDS
         IF(OZONX(IL).GT.OZON(IOX,KL))  THEN
           IOZ(IL) = IOX
           OZON1(IL) = OZON(IOX,KL)
           OZON2(IL) = OZON(IOX+1,KL)
         ENDIF
 160   ENDDO
 170 ENDDO
     DO 190 IL = 1,MERIDS
       X1(IL) = OZONX(IL) - OZON1(IL)
       X2(IL) = OZON2(IL) - OZONX(IL)
       XX = X1(IL) + X2(IL)
       IF(XX.EQ.0.0) THEN
         WRITE(6,*) 'TEST 1 ', MYPE
         WRITE(6,*) IL,KL,X1(IL),X2(IL),OZONX(IL),OZON1(IL),             &
                    OZON2(IL),IOZ(IL)
       ENDIF
       Y1(IL) = 1./(X1(IL) + X2(IL))
       Z1(IL) = COSPHI(IL) - COSP1(IL)
       Z2(IL) = COSP2(IL) - COSPHI(IL)
       ZZ = Z1(IL) + Z2(IL)
       IF(ZZ.EQ.0.0) THEN
         WRITE(6,*) 'TEST 2 ', MYPE
         WRITE(6,*) IL,KL,Z1(IL),Z2(IL),COSPHI(IL),COSP1(IL),            &
                    COSP2(IL),ICOZ(IL)
       ENDIF
       Y2(IL) = 1./ZZ
 190 ENDDO
     DO 200 IL = 1,MERIDS
       X1(IL) = X1(IL)*Y1(IL)
       X2(IL) = X2(IL)*Y1(IL)
       Z1(IL) = Z1(IL)*Y2(IL)
       Z2(IL) = Z2(IL)*Y2(IL)
 200 ENDDO
     DO 220 IR = 1,IRX4
       DO 210 IL = 1,MERIDS
         FOTO(IL,IR) = 0.
 210   ENDDO
 220 ENDDO
     DO 240 IL = 1,MERIDS
       IF(LCOZ(IL))  THEN
         DO 230 IR = 1,IRX4
           P1 = PHOTO(IR,ICOZ(IL)+1,IOZ(IL)+1,KL)*X1(IL) +                 &
                X2(IL)*PHOTO(IR,ICOZ(IL)+1,IOZ(IL),KL)
           P2 = PHOTO(IR,ICOZ(IL),IOZ(IL)+1,KL)*X1(IL) +                   &
             X2(IL)*PHOTO(IR,ICOZ(IL),IOZ(IL),KL)
           ALG =  P1*Z1(IL) + P2*Z2(IL)
           FOTO(IL,IR) = ALG
 230     ENDDO
       ENDIF
 240 ENDDO
!
!  PHOTODISSOCIATION RATES FOR TEMPERATURE DEPENDENT N2O5 CROSS-SECTION 
!  FIT A QUADRATIC TO THE VALUES AT 300K, 250K AND 200K                 
!
     DO 250 IL = 1,MERIDS
       IF(LCOZ(IL))  THEN
         XX = 5.0 - VTEMP(IL)*1000.0
         XX2 = XX*XX
         DO 245 IM = 0,3
           CX1 = FOTO(IL,10+IRX*IM)
           AX1 = (9.*FOTO(IL,8+IRX*IM) - 15.*FOTO(IL,9+IRX*IM) +             &
                  6.*CX1)*0.1
           BX1 = FOTO(IL,9+IRX*IM) - AX1 - CX1
           FOTO(IL,8+IRX*IM) = AX1*XX2 + BX1*XX + CX1
           CX2 = FOTO(IL,21+IRX*IM)      
           AX2 = (9.*FOTO(IL,5+IRX*IM) - 15.*FOTO(IL,20+IRX*IM) +            &
                  6.*CX2)*0.1
           BX2 = FOTO(IL,20+IRX*IM) - AX2 - CX2
           FOTO(IL,5+IRX*IM) = AX2*XX2 + BX2*XX + CX2
           CX3 = FOTO(IL,23+IRX*IM)      
           AX3 = (9.*FOTO(IL,12+IRX*IM) - 15.*FOTO(IL,22+IRX*IM) +           &
                  6.*CX3)*0.1
           BX3 = FOTO(IL,22+IRX*IM) - AX3 - CX3
           FOTO(IL,12+IRX*IM) = AX3*XX2 + BX3*XX + CX3
           CX1 = FOTO(IL,30+IRX*IM)
           AX1 = (9.*FOTO(IL,28+IRX*IM) - 15.*FOTO(IL,29+IRX*IM) +           &
                  6.*CX1)*0.1
           BX1 = FOTO(IL,29+IRX*IM) - AX1 - CX1
           FOTO(IL,28+IRX*IM) = AX1*XX2 + BX1*XX + CX1
 245     ENDDO
       ENDIF
 250 ENDDO
     DO 300 IL = 1,MERIDS
       IF(LCOZ(IL))  THEN  
         DO IM = 0,3
           DO IR = 1,IRX
             IF(FOTO(IL,IR+IM*IRX).GT.10.0) FOTO(IL,IR+IM*IRX) = 10.0
             IF(FOTO(IL,IR+IM*IRX).LT.-100.0) FOTO(IL,IR+IM*IRX) = -100.0
             FOTO(IL,IR+IM*IRX) =EXP(FOTO(IL,IR+IM*IRX))  
           ENDDO
         ENDDO
       ENDIF 
 300 ENDDO
     RETURN
END SUBROUTINE PHOTO_RATES
!/ ----------------------------------------------------------------     
!      
SUBROUTINE ZEN2(ITIME,DT,DRAD,DARG,DELTA,CDX)         
!
!  THIS SUBROUTINE CALCULATES THE COSINE OF THE SOLAR ZENITH ANGLE,     
!  COZEN, ALLOWING FOR LEAP YEARS AND THEELLIPTICITY OF THEEARTH'S    
!  ORBIT.
! ELONG IS THE LONGITUDE IN DEG.E. ,XLAT THE LATITUDE NORTH. GMT IS   
!  THE G.M.T. AND DAY IS THE DAY IN THE MONTH. TIMING INFO. IS GIVEN BY 
!  ITIME.
!
IMPLICIT NONE
INTEGER, intent(in) :: ITIME(6)
REAL, intent(in)    :: DT,DRAD
REAL, intent(inout)   :: DARG,DELTA,CDX

REAL DGMT,DP,D1,D2,DF1,DTY,DR,DET,    &
      DAL1,DAL,CD1,DTD
!
INTEGER  IYEAR,IMON,IDAY
!
    IYEAR = ITIME(1)
    IMON = ITIME(2)
    IDAY = ITIME(3)
    DGMT = ITIME(4) + ITIME(5)/60.0 + (ITIME(6) + 0.5*DT)/3600.0
    IF(DGMT.LT.24.)  GOTO 5
    DGMT = DGMT - 24.
    IDAY = IDAY + 1.
  5 DP = 365.*IYEAR + IDAY + 31.*(IMON-1.)
    IF(IMON - 2)  10,10,20
 10 D1 = (IYEAR - 1.)*0.25
    D2 = (IYEAR - 1.)*1.E-2 + 1.
    DF1 = DP + INT(D1) - INT(0.75*INT(D2))
    GOTO 30
 20 D1 = 4.E-1*IMON + 2.3
    D2 = IYEAR*1.E-2
    DF1 = DP - INT(D1) + INT(25.*D2) - INT(0.75*(INT(D2) + 1.))       
 30 DTY = (DF1 - 693960.)/36525.
    DTD = 100.*DTY
    DR = (2.7969668E02 + 3.6000768925E04*DTY + 3.03E-4*DTY*DTY)*    &
      DRAD
    DET = (-9.32701E01  - 1.42E-1*DTD)*SIN(DR) + (-4.3194497E02 +    &
      3.3E-2*DTD)*COS(DR) + 5.965E02*SIN(2.*DR) - 2.*COS(2.*DR)       &
    + 4.2E00*SIN(3.*DR) + 1.93E01*COS(3.*DR) - 1.28E01*SIN(4.*DR)   
    DAL1 = DET*DRAD/240.0
    DAL = DR - DAL1
!
!  DELTA IS THE SOLAR DECLINATION AND COZEN IS THE COSINE OF THE SOLAR  
!  ZENITH ANGLE
!
    DELTA = ATAN(4.336E-1*SIN(DAL))
    DARG = (DGMT - 12.0)*15.*DRAD + DAL1
    CD1 = DTD*DRAD*360.0
    CDX = 1.00011 + 0.034221*COS(CD1) + 0.001280*SIN(CD1) +            &
          0.000719*COS(2.*CD1) + 0.000077*SIN(2.*CD1)                 
    RETURN
END SUBROUTINE ZEN2
!
SUBROUTINE SEDIMENT(ANAT,AICE,ANOY,AH2O,AHNO3,P_FIELD,LEVELS,     &
                    PRESS,DT,MYPE)      
!
!  CALCULATES SEDIMENTATION RATES OF TYPE I AND TYPE 2 PARTICLES               
!  AND VERTICALLY ADVECTS MODEL NAT AND ICE
!
  INTEGER, intent(in) ::  P_FIELD, LEVELS, MYPE
  REAL, intent(inout) :: ANAT(P_FIELD,LEVELS),AICE(P_FIELD,LEVELS), &
       ANOY(P_FIELD,LEVELS),AH2O(P_FIELD,LEVELS),                   &
       AHNO3(P_FIELD,LEVELS)
  REAL, intent(in) :: PRESS(P_FIELD,LEVELS)
  REAL, intent(in) :: DT   


  REAL  SNATS(P_FIELD,LEVELS),SICES(P_FIELD,LEVELS),                 &
        F1(P_FIELD,LEVELS),F2(P_FIELD,LEVELS), &
        PNAT(P_FIELD),PICE(P_FIELD),PNAT2(P_FIELD),PICE2(P_FIELD),        &
        ANAT2(P_FIELD,LEVELS),AICE2(P_FIELD,LEVELS),                 &
        ANATMAX(P_FIELD),AICEMAX(P_FIELD)
!
!  DATA: V1 IS SEDIMENTATION VELOCITY (M/S)
!  OF ICE PARTICLES V2 IS SEDIMENTATION VELOCITY OF         
!  NAT PARTICLES, ASSUMED RADII R1, R2.
!  AM1, AM2 ARE THE MOLECULAR WEIGHTS AND RHO1, RHO2 ARE THE DENSITIES 
!  OF THE PSCs IN G/CM3
!
!      DATA  V1, V2/1.27E-2, 1.39E-4/
!      DATA  R1, R2/7.0E-6, 0.5E-6/
!      DATA AM1,AM2/18.0,117.0/
!      DATA RHO1,RHO2/0.928, 1.35/
  REAL ::   V1=1.27E-2,  V2=1.39E-4
  REAL ::   R1=7.0E-6,   R2=0.5E-6
  REAL ::  AM1=18.0,    AM2=117.0
  REAL :: RHO1=0.928,  RHO2=1.35
  REAL :: RATIO, TEMP, PFRAC, DZ, CONST, D1, D2, FIXNAT, FIXICE
  INTEGER :: KL, JL

!
!  CALCULATE FRACTION OF NAT PARTICLES USED AS TYPE 2 CORES (F1)               
!  AND FRACTION OF NAT PARTICLES THAT REMAIN AS TYPE 1 CORES (F2)              
!  DETERMINE MAXIMUM NAT AND ICE TO APPLY LIMITERS TO ADVECTED AMOUNTS
!
    ANATMAX(:) = 0.0
    AICEMAX(:) = 0.0
    RATIO = AM1*RHO2/(AM2*RHO1)*(R2/R1)**3
    DO 20 KL = 1,LEVELS 
      DO 10 JL = 1,P_FIELD 
        IF(AICE(JL,KL).LT.1.0E-18) THEN 
          AICE(JL,KL) = 0.       
        ENDIF
        IF(ANAT(JL,KL).LT.1.0E-18)  THEN
          ANAT(JL,KL) = 0.
          F1(JL,KL) = 0.
          F2(JL,KL) = 0.
        ELSE
          F1(JL,KL) = AICE(JL,KL)*RATIO/ANAT(JL,KL)
          IF(F1(JL,KL).GT.1.0)  F1(JL,KL) = 1.0
          F2(JL,KL) = 1.0 - F1(JL,KL)
        ENDIF
        IF(ANAT(JL,KL).GT.ANATMAX(JL)) ANATMAX(JL) = ANAT(JL,KL)
        IF(AICE(JL,KL).GT.AICEMAX(JL)) AICEMAX(JL) = AICE(JL,KL)
 10   ENDDO
 20 ENDDO
!
! VERTICALLY ADVECT NAT AND ICE. NOTE THAT PART OF NAT IS ADVECTED             
! AT TYPE 2 RATE AND THE REMAINDER AT TYPE 1 RATE. CALCULATE DESCENT IN
!  1 TIMESTEP; USE APPROXIMATE VERTICAL DISPLACEMENT BETWEEN LAYERS
!
    TEMP = 195.0
    PNAT(:) = 0.0
    PICE(:) = 0.0
    DO 50 KL = 2,LEVELS
      DO 40 JL = 1,P_FIELD
        PFRAC = PRESS(JL,KL)/PRESS(JL,KL-1)
        DZ = 29.26*TEMP*ALOG(PFRAC)
        CONST = DT/DZ
        D1 = ANAT(JL,KL) - ANAT(JL,KL-1)
        D2 = AICE(JL,KL) - AICE(JL,KL-1)         
        SNATS(JL,KL) = -CONST*D1*(V1*F1(JL,KL) + V2*F2(JL,KL))
        SICES(JL,KL) = -CONST*D2*V1
        PNAT(JL) = PNAT(JL) + PRESS(JL,KL)*ANAT(JL,KL)
        PICE(JL) = PICE(JL) + PRESS(JL,KL)*AICE(JL,KL)
 40   ENDDO!CONTINUE
 50 ENDDO!CONTINUE
!
!  set sedimented nat and ice to zero at top and bottom
!
    SNATS(:,1) = 0.0
    SICES(:,1) = 0.0 
    SNATS(:,LEVELS) = 0.0
    SICES(:,LEVELS) = 0.0
    DO 70 KL = 1,LEVELS
      DO 60 JL = 1,P_FIELD
        ANAT2(JL,KL) = ANAT(JL,KL) + SNATS(JL,KL)
        AICE2(JL,KL) = AICE(JL,KL) + SICES(JL,KL)
!
!  APPLY LIMITERS TO NEW NAT AND ICE
!
        IF(ANAT2(JL,KL).LT.0.0) ANAT2(JL,KL) = 0.0
        IF(AICE2(JL,KL).LT.0.0) AICE2(JL,KL) = 0.0
        IF(ANAT2(JL,KL).GT.ANATMAX(JL)) ANAT2(JL,KL) = ANATMAX(JL)
        IF(AICE2(JL,KL).GT.AICEMAX(JL)) AICE2(JL,KL) = AICEMAX(JL)
 60   ENDDO!CONTINUE
 70 ENDDO!CONTINUE
!
! APPLY MASS FIXER
!
    PNAT2(:) = 0.0
    PICE2(:) = 0.0
    DO 90 KL = 1,LEVELS
      DO 80 JL = 1,P_FIELD
        PNAT2(JL) = PNAT2(JL) + PRESS(JL,KL)*ANAT2(JL,KL)
        PICE2(JL) = PICE2(JL) + PRESS(JL,KL)*AICE2(JL,KL)
 80   ENDDO!CONTINUE
 90 ENDDO!CONTINUE
    DO 110 JL = 1,P_FIELD
      IF(PNAT2(JL).EQ.0.0) THEN 
        FIXNAT = 1.0
      ELSE 
        FIXNAT = PNAT(JL)/PNAT2(JL)
      ENDIF
      IF(PICE2(JL).EQ.0.0) THEN
        FIXICE = 1.0
      ELSE
        FIXICE = PICE(JL)/PICE2(JL)
      ENDIF
      DO 100 KL = 1,LEVELS
        ANAT2(JL,KL) = ANAT2(JL,KL)*FIXNAT 
        AICE2(JL,KL) = AICE2(JL,KL)*FIXICE 
100   ENDDO
110 ENDDO
!
!  ADJUST NOY AND H2O TENDENCY FIELDS
!
    DO KL = 1,LEVELS
      DO JL = 1,P_FIELD
        ANOY(JL,KL) = ANOY(JL,kl) + (ANAT2(JL,KL) - ANAT(JL,KL))/DT 
        AHNO3(JL,KL) = AHNO3(JL,kl) + (ANAT2(JL,KL) - ANAT(JL,KL))/DT 
        AH2O(JL,KL) = AH2O(JL,KL) + (AICE2(JL,KL) - AICE(JL,KL))/DT 
     ENDDO
   ENDDO
   RETURN
END SUBROUTINE SEDIMENT
!
SUBROUTINE GHGS(ITIME,ISC,CH4,AN2O,CH4TR,AN2OTR)
!
!  FINDS THE CONCENTRATIONS OF THE GHGS CH4 AND N2O ACCORDING TO THE TIME 
!  ITIME AND SCENARIO ISC, CHOSEN FROM SRES
!
  INTEGER, intent(in) :: ITIME(6),ISC
  REAL, intent(inout) :: CH4,AN2O,CH4TR,AN2OTR

  REAL TIMEN2O(22),AN2O_A1B(22),AN2O_A2(22),AN2O_B1(22)
  REAL TIMECH4(36),CH4_A1B(36),CH4_A2(36),CH4_B1(36)
  REAL TIME, SECS, X1,X2, Z1, C1, C2
  INTEGER IMON(12), IT
  DATA  IMON/0,31,59,90,120,151,181,212,243,273,304,334/
  DATA TIMEN2O/                                                     &
     1950.5, 1976.5, 1980.5, 1984.5, 1988.5, 1990.5,                &
     1992.5, 1994.5, 1996.5, 1998.5, 2000.5, 2010.5,                &
     2020.5, 2030.5, 2040.5, 2050.5, 2060.5, 2070.5,                &
     2080.5, 2090.5, 2100.5, 2500.5/
  DATA AN2O_A1B/                                                    &
     278.0,  299.0,  301.0,  303.0,  305.0,  309.6,                 &
     311.0,  312.0,  312.4,  314.0,  316.0,  324.0,                 &
     331.0,  338.0,  344.0,  350.0,  356.0,  360.0,                 &
     365.0,  368.0,  372.0,  372.0/
  DATA AN2O_A2/                                                     &
     278.0,  299.0,  301.0,  303.0,  305.0,  309.6,                 &
     311.0,  312.0,  312.4,  314.0,  316.0,  325.0,                 &
     335.0,  347.0,  360.0,  373.0,  387.0,  401.0,                 &
     416.0,  432.0,  447.0,  447.0/
  DATA AN2O_B1/                                                     &
     278.0,  299.0,  301.0,  303.0,  305.0,  309.6,                 &
     311.0,  312.0,  312.4,  314.0,  316.0,  324.0,                 &
     333.0,  341.0,  349.0,  357.0,  363.0,  368.0,                 &
     371.0,  374.0,  375.0,  375.0/
  DATA TIMECH4/                                                     &
     1950.5, 1952.5, 1954.5, 1956.5, 1958.5, 1960.5,                &
     1962.5, 1964.5, 1966.5, 1968.5, 1970.5, 1972.5,                &
     1974.5, 1976.5, 1978.5, 1980.5, 1982.5, 1984.5,                &
     1986.5, 1988.5, 1990.5, 1992.5, 1994.5, 1996.5,                &
     1998.5, 2000.5, 2010.5, 2020.5, 2030.5, 2040.5,                &
     2050.5, 2060.5, 2070.5, 2080.5, 2090.5, 2100.5/
  DATA CH4_A1B/                                                     &
     1147.5, 1163.8, 1182.1, 1202.4, 1224.2, 1247.5,                &
     1272.2, 1298.4, 1326.0, 1355.3, 1386.0, 1417.5,                &
     1449.3, 1481.5, 1514.0, 1547.1, 1580.9, 1614.2,                &
     1644.9, 1671.6, 1694.3, 1714.0, 1721.0, 1728.0,                &
     1745.0, 1760.0, 1871.0, 2026.0, 2202.0, 2337.0,                &
     2400.0, 2386.0, 2301.0, 2191.0, 2078.0, 1974.0/
  DATA CH4_A2/                                                      &
     1147.5, 1163.8, 1182.1, 1202.4, 1224.2, 1247.5,                &
     1272.2, 1298.4, 1326.0, 1355.3, 1386.0, 1417.5,                &
     1449.3, 1481.5, 1514.0, 1547.1, 1580.9, 1614.2,                &
     1644.9, 1671.6, 1694.3, 1714.0, 1721.0, 1728.0,                &
     1745.0, 1760.0, 1861.0, 1997.0, 2163.0, 2357.0,                &
     2562.0, 2779.0, 3011.0, 3252.0, 3493.0, 3731.0/
  DATA CH4_B1/                                                      &
     1147.5, 1163.8, 1182.1, 1202.4, 1224.2, 1247.5,                &
     1272.2, 1298.4, 1326.0, 1355.3, 1386.0, 1417.5,                &
     1449.3, 1481.5, 1514.0, 1547.1, 1580.9, 1614.2,                &
     1644.9, 1671.6, 1694.3, 1714.0, 1721.0, 1728.0,                &
     1745.0, 1760.0, 1827.0, 1891.0, 1927.0, 1919.0,                &
     1881.0, 1836.0, 1797.0, 1741.0, 1663.0, 1574.0/
!
    TIME = ITIME(1) + (IMON(ITIME(2)) + ITIME(3))/365.0
    IF(ISC.EQ.4) TIME = ITIME(1) + 0.5
    IF(TIME.LT.1950.5) TIME = 1950.5
    IF(TIME.GT.2100.5) TIME = 2100.5
    SECS = 365.25*86400.0
    DO IT = 2,22
      IF(TIMEN2O(IT).GT.TIME) THEN
        X1 = TIMEN2O(IT) - TIME
        X2 = TIME - TIMEN2O(IT-1)
        Z1 = X1 + X2
        IF(ISC.EQ.1) THEN
          C1 = AN2O_A1B(IT-1)
          C2 = AN2O_A1B(IT)
        ELSEIF(ISC.EQ.2) THEN
          C1 = AN2O_A2(IT-1)
          C2 = AN2O_A2(IT)
        ELSEIF(ISC.EQ.3) THEN
          C1 = AN2O_B1(IT-1)
          C2 = AN2O_B1(IT)
        ELSEIF(ISC.EQ.4) THEN
          C1 = AN2O_A1B(IT-1)
          C2 = AN2O_A1B(IT)
        ENDIF
        AN2O = (X1*C1 + X2*C2)/Z1
        AN2OTR = (C2 - C1)/(AN2O*SECS*Z1)
        GOTO 10
      ENDIF
    ENDDO
 10 DO IT = 2,36
      IF(TIMECH4(IT).GT.TIME) THEN
        X1 = TIMECH4(IT) - TIME
        X2 = TIME - TIMECH4(IT-1)
        Z1 = X1 + X2
        IF(ISC.EQ.1) THEN
          C1 = CH4_A1B(IT-1)
          C2 = CH4_A1B(IT)
        ELSEIF(ISC.EQ.2) THEN
          C1 = CH4_A2(IT-1)
          C2 = CH4_A2(IT)
        ELSEIF(ISC.EQ.3) THEN
          C1 = CH4_B1(IT-1)
          C2 = CH4_B1(IT)
        ELSEIF(ISC.EQ.4) THEN
          C1 = CH4_A1B(IT-1)
          C2 = CH4_A1B(IT)
        ENDIF
        CH4 = (X1*C1 + X2*C2)/Z1
        CH4TR = (C2 - C1)/(CH4*SECS*Z1)
        GOTO 20
      ENDIF
    ENDDO
 20 AN2O = AN2O*1.0E-9
    CH4 = CH4*1.0E-9
    RETURN
END SUBROUTINE GHGS
!/ ----------------------------------------------------------------            
END MODULE STRAT_CHEM_MOD
