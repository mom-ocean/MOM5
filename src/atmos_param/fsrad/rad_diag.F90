
                     MODULE RAD_DIAG_MOD

!-----------------------------------------------------------------------

!!!   USE   RDPARM_MOD, ONLY:  IMAX, LMAX, LP1, NBLW, NBLY, NBLM
      USE   RDPARM_MOD, ONLY:  LMAX, LP1, NBLW, NBLY, NBLM

      USE   HCONST_MOD, ONLY:  RADCON, RADCON1

      USE LONGWAVE_MOD, ONLY: OSOUR, CSOUR, SS1
      USE LONGWAVE_MOD, ONLY: FLX1E1, GXCTS, FCTSG
      USE LONGWAVE_MOD, ONLY: CLDFAC
      USE LONGWAVE_MOD, ONLY: DELP2, DELP
      USE LONGWAVE_MOD, ONLY: TO3, CO21, EMISS, EMISS2, CTS, EXCTS,  &
                              EXCTSN, E1FLX, CO2SP
      USE LONGWAVE_MOD, ONLY: IBAND, BANDLO, BANDHI

      Use       Fms_Mod, ONLY: write_version_number, mpp_pe, mpp_root_pe, &
                               error_mesg, FATAL


implicit none
private

!-----------------------------------------------------------------------
      character(len=128) :: version = '$Id: rad_diag.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical            :: module_is_initialized = .false.

public RADIAG, RAD_DIAG_init, RAD_DIAG_end


      CONTAINS

!#######################################################################
!#######################################################################

      SUBROUTINE RADIAG  &
          ( PRESS,TEMP,RH2O,RRVCO2,QO3,CAMT,KTOP,KBTM,NCLDS,  &
            HEATRA,GRNFLX,  &
            FSW,DFSW,UFSW,HSW,  &
            KTOPSW,KBTMSW,EMCLD,CUVRF,CIRRF,CIRAB,  &
            SALB,COSZRO,SSOLAR,   ip,jp)

!-----------------------------------------------------------------------
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRESS,TEMP,RH2O
      REAL,    INTENT(IN)                   :: RRVCO2
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: QO3,CAMT
      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)   :: NCLDS

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: HEATRA
      REAL,    INTENT(IN), DIMENSION(:,:)   :: GRNFLX

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD,CUVRF,CIRRF,CIRAB

      REAL,    INTENT(IN), DIMENSION(:,:)   :: SALB,COSZRO,SSOLAR
      INTEGER, INTENT(IN)                 :: ip,jp
!-----------------------------------------------------------------------

      CALL RADPRT (PRESS(ip,jp,:),TEMP(ip,jp,:),RH2O(ip,jp,:),RRVCO2,QO3(ip,jp,:), &
                   CAMT(ip,jp,:),KTOP(ip,jp,:),KBTM(ip,jp,:),NCLDS(ip,jp),         &
                   HEATRA(ip,jp,:),GRNFLX(ip,jp),  &
                   CTS(ip,:),EXCTS(ip,:),    &
                   EMISS(ip,:,:),CLDFAC(ip,:,:),E1FLX(ip,:),        &
                   DELP(ip,:),DELP2(ip,:),EMISS2(ip,:,:),           &
                   SS1(ip,:),TO3(ip,:,:),OSOUR(ip,:),CO21(ip,:,:),  &
                   CSOUR(ip,:),CO2SP(ip,:),     &
                   GXCTS(ip),FLX1E1(ip),        &
                   FCTSG(ip,:),EXCTSN(ip,:,:),  &
                   FSW(ip,jp,:),DFSW(ip,jp,:),UFSW(ip,jp,:),HSW(ip,jp,:),  &
                   KTOPSW(ip,jp,:),KBTMSW(ip,jp,:),EMCLD(ip,jp,:),      &
                   CUVRF(ip,jp,:),CIRRF(ip,jp,:),CIRAB(ip,jp,:),        &
                   SALB(ip,jp),COSZRO(ip,jp),SSOLAR(ip,jp),   ip,jp)

!-----------------------------------------------------------------------

      END SUBROUTINE RADIAG

!#######################################################################
!#######################################################################

      subroutine RAD_DIAG_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      end subroutine RAD_DIAG_init

!#######################################################################
!#######################################################################

      subroutine RAD_DIAG_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

      end subroutine RAD_DIAG_end
!#######################################################################
!#######################################################################

      SUBROUTINE RADPRT (PRESS,TEMP,RH2O,RRVCO2,QO3,  &
                         CAMT,KTOP,KBTM,NCLDS,        &
                         HEATRA,GRNFLX,  &
                         CTS,EXCTS,      &
                         EMISS,CLDFAC,E1FLX,DELP,DELP2,          &
                         EMISS2,SS1,TO3,OSOUR,CO21,CSOUR,CO2SP,  &
                         GXCTS,FLX1E1,       &
                         FCTSG,EXCTSN,       &
                         FSW,DFSW,UFSW,HSW,  &
                         KTOPSW,KBTMSW,EMCLD,CUVRF,CIRRF,CIRAB,  &
                         SALB,COSZRO,SSOLAR,              &
                         IP,JP)

      IMPLICIT NONE
!-----------------------------------------------------------------------
!      Subroutine RADPRT prints out all diagnostic quantities
!      required for diagnosos of radiation quantities.
!-----------------------------------------------------------------------
      REAL,    INTENT(IN), DIMENSION(:) :: PRESS,TEMP,RH2O
      REAL,    INTENT(IN)               :: RRVCO2
      REAL,    INTENT(IN), DIMENSION(:) :: QO3,CAMT
      INTEGER, INTENT(IN), DIMENSION(:) :: KTOP,KBTM
      INTEGER, INTENT(IN)               :: NCLDS

      REAL,    INTENT(IN), DIMENSION(:) :: HEATRA
      REAL,    INTENT(IN)               :: GRNFLX

      REAL,    INTENT(IN), DIMENSION(:)   :: CTS,EXCTS,E1FLX,DELP,DELP2
      REAL,    INTENT(IN), DIMENSION(:,:) :: EMISS,CLDFAC,EMISS2,TO3
      REAL,    INTENT(IN), DIMENSION(:,:) :: CO21
      REAL,    INTENT(IN), DIMENSION(:)   :: SS1,OSOUR,CSOUR,CO2SP

      REAL,    INTENT(IN)                 :: GXCTS,FLX1E1
      REAL,    INTENT(INOUT), DIMENSION(:)   :: FCTSG
      REAL,    INTENT(INOUT), DIMENSION(:,:) :: EXCTSN

      REAL,    INTENT(IN), DIMENSION(:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), DIMENSION(:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:) :: EMCLD,CUVRF,CIRRF,CIRAB

      REAL,    INTENT(IN)                 :: SALB,COSZRO,SSOLAR
      INTEGER, INTENT(IN)                 :: IP,JP
!-----------------------------------------------------------------------
!    ----- DIMENSION FOR LOCAL VARIABLES -----

      REAL  VSUM1(LP1),SUM(LP1),FLXDG(LP1),HTEM(LMAX),HTEM1(LMAX),  &
            HTEM2(LMAX),HTEM3(LMAX),HTEM4(LMAX),FLXNET(LP1),  &
            FTOPN(NBLY),FTOPAC(NBLY),VSUMAC(NBLY),PFLUX(LP1),  &
            CTS1(LMAX),CTS2(LMAX),CTST(LMAX),  &
            OSS(LP1),CSS(LP1),TC(LP1),DTC(LP1),SS2(LP1),  &
            HLWSW(LMAX),FLWSW(LP1),  &
            FSWD(LP1),DFSWD(LP1),UFSWD(LP1)
!-----------------------------------------------------------------------
      INTEGER  K,KP,N,NX,NY,NPRT
      REAL     FTOPC,FTOPE,FTOP,FGRD,FDIFF,SOLARW,QSUM
!-----------------------------------------------------------------------
!****COMPUTE LOCAL VARIABLES TC,DTC,OSS,CSS,SS2

      DO K=1,LP1
         TC(K)=TEMP(K)*TEMP(K)*TEMP(K)*TEMP(K)
      ENDDO

      DO K=2,LP1
         DTC(K)=TC(K)-TC(K-1)
         SS2(K)=SS1(K)-SS1(K-1)
         CSS(K)=CSOUR(K)-CSOUR(K-1)
         OSS(K)=OSOUR(K)-OSOUR(K-1)
      ENDDO

!***HTEM1 = EMISSIVITY HEATING RATE FOR 0-160,1200-2200 CM-1 BAND

      DO K=1,LP1
         SUM(K)=0.
         DO KP=1,LMAX
            VSUM1(KP)=DTC(KP+1)*EMISS(KP+1,K)*CLDFAC(KP+1,K)
         ENDDO
         DO KP=1,LMAX
            SUM(K)=SUM(K)+VSUM1(KP)
         ENDDO
      ENDDO

      DO K=1,LP1
         FLXDG(K)=SUM(K)+TC(1)*E1FLX(K)*CLDFAC(K,1)
      ENDDO

      DO K=1,LMAX
         HTEM1(K)=RADCON*(FLXDG(K+1)-FLXDG(K))*DELP(K)
      ENDDO

!***HTEM2 = EMISSIVITY HEATING RATE FOR 800-990,1070-1200 CM-1 BAND

      DO K=1,LP1
         SUM(K)=0.
         DO KP=1,LMAX
            VSUM1(KP)=SS2(KP+1)*EMISS2(KP+1,K)*CLDFAC(KP+1,K)
         ENDDO
         DO KP=1,LMAX
            SUM(K)=SUM(K)+VSUM1(KP)
         ENDDO
      ENDDO

      DO K=1,LP1
         FLXDG(K)=SUM(K)+SS1(1)*EMISS2(K,1)*CLDFAC(K,1)
      ENDDO

      DO K=1,LMAX
         HTEM2(K)=RADCON*(FLXDG(K+1)-FLXDG(K))*DELP(K)
      ENDDO

!***HTEM3 = EMISSIVITY HEATING RATE FOR 990-1070 CM-1 BAND

      DO K=1,LP1
         SUM(K)=0.
         DO KP=1,LMAX
            VSUM1(KP)=OSS(KP+1)*TO3(KP+1,K)*CLDFAC(KP+1,K)
         ENDDO
         DO KP=1,LMAX
            SUM(K)=SUM(K)+VSUM1(KP)
         ENDDO
      ENDDO

      DO K=1,LP1
         FLXDG(K)=SUM(K)+OSOUR(1)*TO3(K,1)*CLDFAC(K,1)
      ENDDO

      DO K=1,LMAX
         HTEM3(K)=RADCON*(FLXDG(K+1)-FLXDG(K))*DELP(K)
      ENDDO

!***HTEM4 = EMISSIVITY HEATING RATE FOR 560-800 CM-1 BAND

      DO K=1,LP1
         SUM(K)=0.
         DO KP=1,LMAX
            VSUM1(KP)=CSS(KP+1)*CO21(KP+1,K)*CLDFAC(KP+1,K)
         ENDDO
         DO KP=1,LMAX
            SUM(K)=SUM(K)+VSUM1(KP)
         ENDDO
      ENDDO

      DO K=1,LP1
         FLXDG(K)=SUM(K)+CSOUR(1)*CO2SP(K)*CLDFAC(K,1)
      ENDDO

      DO K=1,LMAX
         HTEM4(K)=RADCON*(FLXDG(K+1)-FLXDG(K))*DELP(K)
      ENDDO

!***HTEM = TOTAL APPROXIMATE HEATING RATE  (Q (APPROX))

      DO K=1,LMAX
         HTEM(K)=HTEM1(K)+HTEM2(K)+HTEM3(K)+HTEM4(K)
      ENDDO

!***COMPUTE FLUX PRESSURE LEVELS

      PFLUX(1)=0.
      DO K=2,LMAX
!!cgs   IF (PRESS(K-1).LT.100.) THEN
        IF (PRESS(K-1).LT.10.) THEN
           PFLUX(K)=2.0*PRESS(K-1)-PFLUX(K-1)
        ELSE
           PFLUX(K)=0.50*(PRESS(K)+PRESS(K-1))
        END IF
      ENDDO
      PFLUX(LP1)=PRESS(LP1)

!***COMPUTE FLUXES AT TOP, GROUND AND NET FLUX AT ALL LEVELS

      FTOPC=GXCTS*1.E-3
      FTOPE=FLX1E1*1.E-3
      FTOP=FTOPC+FTOPE
      FGRD=GRNFLX*1.E-3
      FDIFF=FTOP-FGRD
      FLXNET(1)=FTOP

      DO K=2,LP1
         FLXNET(K)=FLXNET(K-1)+HEATRA(K-1)*DELP2(K-1)*RADCON1*1.E-3
      ENDDO

      DO K=1,LP1
         FSWD(K)=FSW(K)*1.E-3
         DFSWD(K)=DFSW(K)*1.E-3
         UFSWD(K)=UFSW(K)*1.E-3
      ENDDO

!***COMPUTE NET RADIATIVE HEATING AND NET FLUX (UP-DOWN)

      SOLARW=SSOLAR*6.97667E5*1.E-3
      DO K=1,LP1
         FLWSW(K)=FLXNET(K)+FSWD(K)
      ENDDO
      DO K=1,LMAX
         HLWSW(K)=HSW(K)+HEATRA(K)
      ENDDO

!***THE CODE BELOW IS FOR DIAGNOSIS OF CTS QUANTITIES***

      DO N=1,NBLM
         QSUM=0.
         DO K=1,LMAX
            QSUM=QSUM+EXCTSN(K,N)/(DELP(K)*RADCON)
         ENDDO
         FTOPN(N)=FCTSG(N)-QSUM
      ENDDO

!***THIS STMT. ACCOUNTS FOR 4.3 UM BAND CONTRIB.

         FTOPN(NBLY)=0.
      DO N=1,NBLY
         FTOPN(N)=FTOPN(N)*1.E-3
      ENDDO

         FTOPAC(1)=FTOPN(1)
      DO N=2,NBLY
         FTOPAC(N)=FTOPAC(N-1)+FTOPN(N)
      ENDDO

!***THIS STMT. ACCOUNTS FOR 4.3 UM BAND CONTRIB.
      FCTSG(NBLY)=0.
!***THIS LOOP SETS EXCTS CONTRIB. OF BAND 15 (4.3UM) TO ZERO

      DO K=1,LMAX
         EXCTSN(K,NBLY)=0.
      ENDDO

      DO N=1,NBLY
         FCTSG(N)=FCTSG(N)*1.E-3
      ENDDO

         VSUMAC(1)=FCTSG(1)
      DO N=2,NBLY
         VSUMAC(N)=VSUMAC(N-1)+FCTSG(N)
      ENDDO

!***APPROXIMATE CTS HEATING RATES***

      DO K=1,LMAX
         CTS1(K)=RADCON*DELP(K)*(CSOUR(K)*  &
            (CO2SP(K+1)*CLDFAC(K+1,1)-CO2SP(K)*CLDFAC(K,1)))
         CTS2(K)=RADCON*DELP(K)*(OSOUR(K)*  &
          (TO3(K+1,1)*CLDFAC(K+1,1)-TO3(K,1)*CLDFAC(K,1)))
         CTST(K)=CTS1(K)+CTS2(K)+CTS(K)
      ENDDO

!-----------------------------------------------------------------------
!****PRINT STATEMENTS FOLLOW***
!***PRINT LAT. POINT AND CLOUD DATA (IF ANY) AND SFC ALBEDO

      PRINT 5650,IP,JP
 5650 FORMAT (//,'  RADIATION RESULTS FOR IP,JP= ',2I5)

      PRINT 5651,NCLDS
 5651 FORMAT (/,' NO. CLOUDS= ',I2)

      IF (NCLDS .ne. 0) THEN
         PRINT 5652
 5652    FORMAT (37X,' LW CLOUD DATA'/,' CLD. NO',8X,'CLD. AMT.',2X,  &
                 'CLD. EMISS',4X,'CLD TOP INDEX',2X,'CLD BOT INDEX')
         PRINT 5653,(N,CAMT(N+1),EMCLD(N+1),KTOP(N+1),  &
                     KBTM(N+1),N=1,NCLDS)
 5653    FORMAT (I5,7X,F12.6,F12.6,I10,I15)
         PRINT 6652
 6652    FORMAT (37X,' SW CLOUD DATA'/,' CLD. NO',8X,'CLD. AMT.',2X,  &
               'CLD TOP INDEX',2X,'CLD BOT INDEX',2X,'VIS. REFL',3X,  &
               ' IR REFL',4X,' IR ABS.')
         PRINT 6653,(N,CAMT(N+1),KTOPSW(N+1),KBTMSW(N+1),  &
                     CUVRF(N+1),CIRRF(N+1),CIRAB(N+1),N=1,NCLDS)
 6653    FORMAT (I5,7X,F12.6,I8,I15,6X,3F12.6)
      ENDIF

!!    NN=NCLDS+2
!!    PRINT 7653, CUVRF(NN),CIRRF(NN)
      PRINT 7653, SALB,SALB
 7653 FORMAT (/,10X,'VIS. SFC. ALBEDO=',F12.6,' IR SFC. ALBEDO=',F12.6)

!***PRINT CO2 AMOUNT***

      PRINT 5654, RRVCO2
 5654 FORMAT (/,' CO2 VOL. MIXING RATIO= ',F14.6,/)

!***PRINT SOLAR INPUT,ZENITH ANGLE***

      PRINT 6654, SOLARW,COSZRO
 6654 FORMAT (/,' INCOMING SOLAR FLUX =',F12.6,' W/M**2',/,  &
                ' COS(AZIMUTH)=',F12.6)

!***PRINT SOLAR INPUT,ZENITH ANGLE,DAY FRACTION***
!     PRINT 6654, SOLARW,COSZRO,TAUDAR
!6654 FORMAT (/,' INCOMING SOLAR FLUX =',F12.6,' W/M**2',/,  &
!               ' COS(AZIMUTH)=',F12.6,10X,' FRACTION SUNUP=',F12.6)

!***PRINT INPUT DATA AND OVERALL HEATING RATES AND FLUXES

      PRINT 5755
 5755 FORMAT (/,20X,' LW HEATING RATES AND FLUXES',/)

      PRINT 5655
 5655 FORMAT ('  LVL',' PRESSURE   ',4X,' TEMP.     ','H2O MMR',5X,  &
              'O3 MMR',7X,'HEAT RATE',2X,'NET FLUX',3X,'FLUX PRESS.')

      PRINT 5555, (K,PRESS(K),TEMP(K),RH2O(K),QO3(K),  &
                     HEATRA(K),FLXNET(K),PFLUX(K),K=1,LMAX)
      PRINT 5556, PRESS(LP1),TEMP(LP1),FLXNET(LP1),PFLUX(LP1)
 5555 FORMAT (I4,E13.6,F12.4,2E12.5,2F12.6,E13.6)
 5556 FORMAT (4X,E13.6,F12.4,36X,F12.6,E13.6)

      PRINT 6755
      PRINT 6655
      PRINT 6555, (K,PRESS(K),HSW(K),FSWD(K),DFSWD(K),  &
                     UFSWD(K),PFLUX(K),K=1,LMAX)
      PRINT 6556, PRESS(LP1),FSWD(LP1),DFSWD(LP1),  &
                  UFSWD(LP1),PFLUX(LP1)
 6755 FORMAT (/,20X,' SW HEATING RATES AND FLUXES',/)
 6655 FORMAT ('  LVL',' PRESSURE    ',3X,'HEAT RATE',2X,'NET FLUX',  &
                4X,'DN FLUX',6X,'UP FLUX',3X,'FLUX PRESS.')
 6556 FORMAT (4X,E13.6,12X,3F12.6,E13.6)
 6555 FORMAT (I4,E13.6,4F12.6,E13.6)

      PRINT 7755
      PRINT 7655
      PRINT 7555, (K,PRESS(K),HLWSW(K),FLWSW(K),PFLUX(K),K=1,LMAX)
      PRINT 7556,    PRESS(LP1),FLWSW(LP1),PFLUX(LP1)
 7755 FORMAT (/,20X,' COMBINED HEATING RATES AND FLUXES',/)
 7655 FORMAT ('  LVL',' PRESSURE    ',4X,'HEAT RATE',2X,'NET FLUX',  &
               3X,'FLUX PRESS.')
 7556 FORMAT (4X,E13.6,12X,F12.6,E13.6)
 7555 FORMAT (I4,E13.6,2F12.6,E13.6)

!***PRINT APPROXIMATE HEATING RATES

      PRINT 5659
 5659 FORMAT (/,37X,'APPROXIMATE HEATING RATES        (Q(APPROX))'/  &
                '  LVL',' PRESSURE   ',5X,'  0-160,1200-2200 ',      &
                ' 800-990,1070-1200','     990-1070     ',           &
                '      560-800     ','       TOTAL')

      PRINT 5559, (K,PRESS(K),HTEM1(K),HTEM2(K),HTEM3(K),HTEM4(K),  &
                     HTEM(K),K=1,LMAX)
 5559 FORMAT (I4,E13.6,5F18.6)

      PRINT 5561
 5561 FORMAT (/,37X,'APPROXIMATE CTS HEATING RATES',  &
              /,'  LVL',' PRESSURE',7X,' H2O BANDS    ',  &
                ' 15 UM BAND   ',' 9.6 UM BAND  ',' TOTAL')

      PRINT 5562, (K,PRESS(K),CTS(K),CTS1(K),CTS2(K),CTST(K),K=1,LMAX)
 5562 FORMAT (I4,E13.6,4F14.6)

      PRINT 5570
 5570 FORMAT (/,37X,'EXACT CTS HEATING RATES, BY BAND',  &
              /,'  LVL',' PRESSURE   ','    TOTAL    ',  &
              5X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7',/)

      PRINT 5573, (K,PRESS(K),EXCTS(K),EXCTSN(K,1),EXCTSN(K,2),  &
                     EXCTSN(K,3),EXCTSN(K,4),EXCTSN(K,5),  &
                     EXCTSN(K,6),EXCTSN(K,7),K=1,LMAX)
 5573 FORMAT (I4,E13.6,8F12.6)

      PRINT 5572
 5572 FORMAT ('  LVL',' PRESSURE   ',7X,'8',11X,'9',10X,'10',10X,  &
                '11',10X,'12',10X,'13',10X,'14',10X,'15')

      PRINT 5573, (K,PRESS(K),EXCTSN(K,8),EXCTSN(K,9),EXCTSN(K,10),  &
                     EXCTSN(K,11),EXCTSN(K,12),EXCTSN(K,13),  &
                     EXCTSN(K,14),EXCTSN(K,15),K=1,LMAX)

!***PRINT FLUXES

      PRINT 5558, FTOPC,FTOPE,FTOP,FGRD,FDIFF
 5558 FORMAT ( 40X,'   FLUXES',  &
             /,' FLUX AT TOP,160-1200 CM-1       =',F14.6,' W/M**2',  &
             /,' FLUX AT TOP,0-160,1200-2200 CM-1=',F14.6,' W/M**2',  &
             /,' FLUX AT TOP,0-2200 CM-1         =',F14.6,' W/M**2',  &
             /,' NET FLUX AT GROUND,0-2200 CM-1  =',F14.6,' W/M**2',  &
             /,' NET FLUX DIFFERENCE,0-2200 CM-1 =',F14.6,' W/M**2')

      PRINT 5633
 5633 FORMAT (40X,'CTS FLUXES',  &
             /,1X,'BAND NO',8X,'LOFREQ',9X,'HIFREQ',9X,'F(1)',11X,  &
                  'ACCUM. F(1)',4X,'CTS F(GRD)',5X,'ACCUM. CTS F(GRD)')

!***PRINTOUT FOR 8 COMBINED BANDS,160-560 CM-1

      DO 767 NY=1,8
      NPRT=1
        DO 768 NX=1,40
         IF (IBAND(NX).EQ.NY) THEN
           IF (NPRT.EQ.1) THEN
              PRINT 5644,NY,BANDLO(NX+16),BANDHI(NX+16),  &
                         FTOPN(NY),FTOPAC(NY),FCTSG(NY),VSUMAC(NY)
              NPRT=0
           ELSE
              PRINT 5646,BANDLO(NX+16),BANDHI(NX+16)
           ENDIF
         ENDIF
768   CONTINUE
767   CONTINUE

!***PRINTOUT FOR REMAINING BANDS

      DO 769 NY=9,NBLM
        PRINT 5644, NY,BANDLO(NY+48),BANDHI(NY+48),FTOPN(NY),  &
                       FTOPAC(NY),FCTSG(NY),VSUMAC(NY)
769   CONTINUE
      NY=NBLY
        PRINT 5644, NY,BANDLO(NBLW),BANDHI(NBLW),FTOPN(NY),  &
                       FTOPAC(NY),FCTSG(NY),VSUMAC(NY)
5644  FORMAT (I11,6F15.6)
5646  FORMAT (11X,2F15.6)

!-----------------------------------------------------------------------

      END SUBROUTINE RADPRT

!#######################################################################
!#######################################################################

                 END MODULE RAD_DIAG_MOD

