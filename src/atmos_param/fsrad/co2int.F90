
      Module co2int_mod

!-----------------------------------------------------------------------
!
!       CO2INT INTERPOLATES CARBON DIOXIDE TRANSMISSION FUNCTIONS
!  FROM THE 109 LEVEL GRID,FOR WHICH THE TRANSMISSION FUNCTIONS
!  HAVE BEEN PRE-CALCULATED, TO THE GRID STRUCTURE SPECIFIED BY THE
!  USER.
!
!        METHOD: 
!
!      CO2INT IS EMPLOYABLE FOR TWO PURPOSES: 1) TO OBTAIN TRANSMIS-
!  SIVITIES BETWEEN ANY 2 OF AN ARRAY OF USER-DEFINED PRESSURES; AND
!  2) TO OBTAIN LAYER-MEAN TRANSMISSIVITIES BETWEEN ANY 2 OF AN ARRAY
!  OF USER-DEFINED PRESSURE LAYERS.TO CLARIFY THESE TWO PURPOSES,SEE
!  THE DIAGRAM AND DISCUSSION BELOW.
!      CO2INT MAY BE USED TO EXECUTE ONLY ONE PURPOSE AT ONE TIME.
!
!     LET P BE AN ARRAY OF USER-DEFINED PRESSURES
!     AND PD BE USER-DEFINED PRESSURE LAYERS.
!
!       - - - - - - - - -   PD(I-1) ---
!                                     !
!       -----------------   P(I)      !  PRESSURE LAYER I  (PLM(I))
!                                     !
!       - - - - - - - - -   PD(I)  ---
!                                     !
!       -----------------   P(I+1)    !  PRESSURE LAYER I+1 (PLM(I+1))
!                                     !
!       - - - - - - - - -   PD(I+1)---
!            ...                          (THE NOTATION USED IS
!            ...                          CONSISTENT WITH THE CODE)
!            ...
!      - - - - - - - - -    PD(J-1)
!
!      -----------------    P(J)
!
!      - - - - - - - - -    PD(J)
!
!      PURPOSE 1:   THE TRANSMISSIVITY BETWEEN SPECIFIC PRESSURES
!      P(I) AND P(J) ,TAU(P(I),P(J))  IS COMPUTED BY THIS PROGRAM.
!      IN THIS MODE,THERE IS NO REFERENCE TO LAYER PRESSURES PD
!      (PD,PLM ARE NOT INPUTTED).
!
!      PURPOSE 2:   THE LAYER-MEAN TRANSMISSIVITY BETWEEN A LAYER-
!      MEAN PRESSURE PLM(J) AND PRESSURE LAYER I IS GIVEN BY
!         TAULM(PLM(I),PLM(J)). IT IS COMPUTED BY THE INTEGRAL
!
!                           PD(I)
!                           ----
!             1             !
!        -------------  *   !   TAU ( P',PLM(J) )  DP'
!        PD(I)-PD(I-1)      !
!                        ----
!                        PD(I-1)
!
!           THE LAYER-MEAN PRESSURE PLM(I) IS SPECIFIED BY THE USER.
!        FOR MANY PURPOSES,PLM WILL BE CHOSEN TO BE THE AVERAGE
!        PRESSURE IN THE LAYER-IE,PLM(I)=0.5*(PD(I-1)+PD(I)).
!           FOR LAYER-MEAN TRANSMISSIVITIES,THE USER THUS INPUTS
!        A PRESSURE ARRAY (PD) DEFINING THE PRESSURE LAYERS AND AN
!        ARRAY (PLM) DEFINING THE LAYER-MEAN PRESSURES.THE CALCULATION
!        DOES NOT DEPEND ON THE P ARRAY USED FOR PURPOSE 1 (P IS NOT
!        INPUTTED).
!
!            THE FOLLOWING PARAGRAPHS DEPICT THE UTILIZATION OF THIS
!       CODE WHEN USED TO COMPUTE TRANSMISSIVITIES BETWEEN SPECIFIC
!       PRESSURES. LATER PARAGRAPHS DESCRIBE ADDITIONAL FEATURES NEEDED
!       FOR LAYER-MEAN TRANSMISSIVITIES.
!
!          FOR A GIVEN CO2 MIXING RATIO AND STANDARD TEMPERATURE
!      PROFILE,A TABLE OF TRANSMISSION FUNCTIONS FOR A FIXED GRID
!     OF ATMOSPHERIC PRESSURES HAS BEEN PRE-CALCULATED.
!      THE STANDARD TEMPERATURE PROFILE IS COMPUTED FROM THE US
!     STANDARD ATMOSPHERE (1977) TABLE.ADDITIONALLY, THE
!     SAME TRANSMISSION FUNCTIONS HAVE BEEN PRE-CALCULATED FOR A
!     TEMPERATURE PROFILE INCREASED AND DECREASED (AT ALL LEVELS)
!     BY 25 DEGREES.
!         THIS PROGRAM READS IN THE PRESPECIFIED TRANSMISSION FUNCTIONS
!     AND A USER-SUPPLIED PRESSURE GRID (P(I)) AND CALCULATES TRANS-
!     MISSION FUNCTIONS ,TAU(P(I),P(J)), FOR ALL P(I)'S AND P(J)'S.
!     A LOGARITHMIC INTERPOLATION SCHEME IS USED.
!         THIS METHOD IS REPEATED FOR THE THREE TEMPERATURE PROFILES
!     GIVEN ABOVE .THEREFORE OUTPUTS FROM THE PROGRAM ARE THREE TABLES
!     OF TRANSMISSION FUNCTIONS FOR THE USER-SUPPLIED PRESSURE GRID.
!     THE EXISTENCE OF THE THREE TABLES PERMITS SUBSEQUENT INTERPO-
!     LATION TO A USER-SUPPLIED TEMPERATURE PROFILE USING THE METHOD
!     DESCRIBED IN THE REFERENCE.SEE LIMITATIONS SECTION IF THE
!     USER DESIRES TO OBTAIN ONLY 1 TABLE OF TRANSMISSIVITIES.
!
!     MODIFICATIONS FOR LAYER-MEAN TRANSMISSIVITIES: 
!          THE PRESSURES INPUTTED ARE THE LAYER-MEAN PRESSURES,PD,
!     AND THE LAYER-MEAN PRESSURES ,PLM. A SERIES OF TRANSMISSIVITIES
!     (TAU(P',PLM(J)) ARE COMPUTED AND THE INTEGRAL GIVEN IN THE
!     DISCUSSION OF PURPOSE 2 IS COMPUTED.FOR PLM(I) NOT EQUAL TO
!     PLM(J) SIMPSON'S RULE IS USED WITH 5 POINTS. IF PLM(I)=PLM(J)
!     (THE "NEARBY LAYER" CASE) A 49-POINT QUADRATURE IS USED FOR
!     GREATER ACCURACY.THE OUTPUT IS IN TAULM(PLM(I),PLM(J)).
!        NOTE: 
!     TAULM IS NOT A SYMMETRICAL MATRIX. FOR THE ARRAY ELEMENT
!     TAULM(PLM(I),PLM(J)),THE INNER(FIRST,MOST RAPIDLY VARYING)
!     DIMENSION IS THE VARYING LAYER-MEAN PRESSURE,PLM(I);THE OUTER
!     (SECOND) DIMENSION IS THE FIXED LAYER-MEAN PRESSURE PLM(J).
!     THUS THE ELEMENT TAULM(2,3) IS THE TRANSMISSION FUNCTION BETWEEN
!     THE FIXED PRESSURE PLM(3) AND THE PRESSURE LAYER HAVING AN AVERAGE
!     PRESSURE OF PLM(2).
!         ALSO NOTE THAT NO QUADRATURE IS PERFORMED OVER THE LAYER
!     BETWEEN THE SMALLEST NONZERO PRESSURE AND ZERO PRESSURE;
!     TAULM IS TAULM(0,PLM(J)) IN THIS CASE,AND TAULM(0,0)=1.
!
!
!             REFERENCE: 
!         S.B.FELS AND M.D.SCHWARZKOPF,"AN EFFICIENT,ACCURATE
!     ALGORITHM FOR CALCULATING CO2 15 UM BAND COOLING RATES",JOURNAL
!     OF GEOPHYSICAL RESEARCH,VOL.86,NO. C2, PP.1205-1232,1981.
!        MODIFICATIONS TO THE ALGORITHM HAVE BEEN MADE BY THE AUTHORS;
!     CONTACT S.B.F.OR M.D.S. FOR FURTHER DETAILS.A NOTE TO J.G.R.
!     IS PLANNED TO DOCUMENT THESE CHANGES.
!
!            AUTHOR:    M.DANIEL SCHWARZKOPF
!
!            DATE:      14 JULY 1983
!
!            ADDRESS: 
!
!                      G.F.D.L.
!                      P.O.BOX 308
!                      PRINCETON,N.J.08540
!                      U.S.A.
!            TELEPHONE:  (609) 452-6521
!
!            INFORMATION ON TAPE: THIS SOURCE IS THE FIRST FILE
!        ON THIS TAPE.THE SIX FILES THAT FOLLOW ARE CO2 TRANS-
!        MISSIVITIES FOR THE 500-850 CM-1 INTERVAL FOR CO2
!        CONCENTRATIONS OF 330 PPMV (1X) ,660 PPMV (2X), AND
!        1320 PPMV (4X). THE FILES ARE ARRANGED AS FOLLOWS: 
!          FILE 2   1X,CONSOLIDATED USING B(250) WEIGHTING FCTN.
!          FILE 3   1X,CONSOLIDATED WITH NO WEIGHTING FCTN.
!          FILE 4   2X,CONSOLIDATED USING B(250) WEIGHTING FCTN.
!          FILE 5   2X,CONSOLIDATED WITH NO WEIGHTING FCTN.
!          FILE 6   4X,CONSOLIDATED USING B(250) WEIGHTING FCTN.
!          FILE 7   4X,CONSOLIDATED WITH NO WEIGHTING FCTN.
!            FILES 2,4,6 ARE RECOMMENDED FOR USE IN OBTAINING
!        TRANSMISSION FUNCTIONS FOR USE IN HEATING RATE
!        COMPUTATIONS;THEY CORRESPOND TO THE TRANSMISSIVITIES
!        DISCUSSED IN THE 1980 PAPER.FILES 3,5,7 ARE PROVIDED
!        TO FACILITATE COMPARISON WITH OBSERVATION AND WITH OTHER
!        CALCULATIONS.
!
!            PROGRAM LANGUAGE: FORTRAN 1977,INCLUDING PARAMETER
!        AND PROGRAM STATEMENTS.THE PROGRAM IS WRITTEN ON A
!        CYBER 170-730.SEE THE SECTION ON LIMITATIONS FOR
!        ADAPTATIONS TO OTHER MACHINES.
!
!
!            PARAMETER INPUTS: 
!     A) NLEVLS    : NLEVLS IS AN (INTEGER) PARAMETER DENOTING
!        THE NUMBER OF NONZERO PRESSURE LEVELS FOR PURPOSE 1
!        OR THE NUMBER OF NONZERO LAYER PRESSURES NEEDED TO
!        SPECIFY THE PRESSURE LAYERS(PURPOSE 2) IN THE OUTPUT
!        GRID. FOR EXAMPLE,IN PURPOSE 1,IF P=0,100,1000,NLEVLS=2.
!        IF,IN PURPOSE 2,PD=0,100,500,1000,THE NUMBER OF NONZERO
!        PRESSURE LAYERS=2,SO NLEVLS=2
!           IN THE CODE AS WRITTEN,NLEVLS=40; THE USER SHOULD
!        CHANGE THIS VALUE TO A USER-SPECIFIED VALUE.
!     B) NLP1,NLP2 : INTEGER PARAMETERS DEFINED AS: NLP1=NLEVLS+1;
!        NLP2=NLEVLS+2.
!           SEE LIMITATIONS FOR CODE MODIFICATIONS IF PARAMETER
!        STATEMENTS ARE NOT ALLOWED ON YOUR MACHINE.
!
!            INPUTS: 
!
!     A) TRANSA    : THE 109X109 GRID OF TRANSMISSION FUNCTIONS
!            TRANSA IS A  REAL ARRAY.
!
!           TRANSA  IS READ FROM FILE 20. THIS FILE CONTAINS 3
!     RECORDS,AS FOLLOWS: 
!        1)   TRANSA, STANDARD TEMPERATURE PROFILE
!        3)   TRANSA, STANDARD TEMPERATURES + 25 DEG
!        5)   TRANSA, STANDARD TEMPERATURES - 25 DEG
!
!     B)   NMETHD: AN INTEGER WHOSE VALUE IS EITHER 1 (IF CO2INT IS
!       TO BE USED FOR PURPOSE 1) OR 2 (IF CO2INT IS TO BE USED FOR
!       PURPOSE 2).
!
!     C)     P,PD,PLM : 
!          P IS A REAL ARRAY (LENGTH NLP1) SPECIFYING THE PRESSURE
!       GRID AT WHICH TRANSMISSION FUNCTIONS ARE TO BE COMPUTED FOR
!       PURPOSE 1.THE DIMENSION  OF P IS  IN MILLIBARS.THE
!       FOLLOWING LIMITATIONS WILL BE EXPLAINED MORE
!       IN THE SECTION ON LIMITATIONS: P(1) MUST BE ZERO; P(NLP1),THE
!       LARGEST PRESSURE, MUST NOT EXCEED 1165 MILLIBARS.
!         PD IS A REAL ARRAY (LENGTH NLP2) SPECIFYING THE PRESSURE
!       LAYERS FOR WHICH LAYER-AVERAGED TRANSMISSION FUNCTIONS ARE
!       TO BE COMPUTED.THE DIMENSION OF PD IS MILLIBARS.THE LIMITATIONS
!       FOR PD ARE THE SAME AS FOR P,AND ARE GIVEN IN THE SECTION ON
!       LIMITATIONS.
!         PLM IS A REAL ARRAY (LENGTH NLP2) SPECIFYING THE LAYER-MEAN
!       PRESSURES. THE DIMENSION OF PLM IS MILLIBARS. THE LIMITATIONS
!       FOR PLM ARE THE SAME AS FOR P,AND ARE GIVEN IN THE SECTION ON
!       LIMITATIONS.PD IS READ IN BEFORE PLM.
!
!          NOTE: AGAIN,WE NOTE THAT THE USER WILL INPUT EITHER P (FOR
!       PURPOSE 1) OR PD AND PLM(FOR PURPOSE 2) BUT NOT BOTH.
!
!
!
!
!           LIMITATIONS: 
!     1)       P(1)=0.,PD(1)=0.,PLM(1)=0. THE TOP PRESSURE LEVEL
!       MUST BE ZERO,OR THE TOP PRESSURE LAYER MUST BE BOUNDED BY ZERO.
!       THE TOP LAYER-MEAN PRESSURE (PLM(1)) MUST BE ZERO; NO
!       QUADRATURE IS DONE ON THE TOP PRESSURE LAYER.EVEN IF ONE IS
!       NOT INTERESTED IN THE TRANSMISSION FUNCTION BETWEEN 0 AND P(J),
!       ONE MUST INCLUDE SUCH A LEVEL.
!     2)      PD(NLP2)=P(NLP1) IS LESS THAN OR EQUAL TO 1165 MB.
!       EXTRAPOLATION TO HIGHER PRESSURES IS NOT POSSIBLE.
!     3)      IF PROGRAM IS NOT PERMITTED ON YOUR COMPILER,
!       SIMPLY DELETE THE LINE.
!     4)      IF PARAMETER IS NOT PERMITTED,DO THE FOLLOWING: 
!            1) DELETE ALL PARAMETER STATEMENTS IN CO2INT
!            2) AT THE POINT WHERE NMETHOD IS READ IN,ADD: 
!                READ (5,202) NLEVLS
!                NLP1=NLEVLS+1
!                NLP2=NLEVLS+2
!            3) CHANGE DIMENSION AND/OR COMMON STATEMENTS DEFINING
!              ARRAYS TRNS,DELTA,P,PD,TRNFCT,PS,PDS,PLM IN CO2INT.
!              THE NUMERICAL VALUE OF (NLEVLS+1) SHOULD BE INSERTED
!              IN DIMENSION OR COMMON STATEMENTS FOR TRNS,DELTA,
!              P,TRNFCT,PS,PLM; THE NUMERICAL VALUE OF (NLEVLS+2)
!              IN DIMENSION OR COMMON STATEMENTS FOR PD,PDS.
!      5)    PARAMETER (NLEVLS=40) AND THE OTHER PARAMETER
!       STATEMENTS ARE WRITTEN IN CDC FORTRAN; ON OTHER MACHINES THE
!       SAME STATEMENT MAY BE WRITTEN DIFFERENTLY,FOR EXAMPLE AS
!       PARAMETER   NLEVLS=40
!      6) "REAL(KIND=8)" IS USED INSTEAD OF "DOUBLE PRECISION" OR 
!           "REAL*8" FOR CODE PORTABILITY.
!       REQUIREMENTS OF CDC FORTAN.
!      7) THE STATEMENT "DO 400 KKK=1,3" CONTROLS THE NUMBER OF
!       TRANSMISSIVITY OUTPUT MATRICES PORDUCED BY THE PROGRAM.TO
!       PRODUCE 1 OUTPUT MATRIX,DELETE THIS STATEMENT.
!
!-----------------------------------------------------------------------

      Use fs_profile_mod, ONLY:  pd1013,plm1013,pd810,plm810
      Use        fms_mod, ONLY:  ERROR_MESG, FATAL, WARNING, &
                                 mpp_pe, mpp_root_pe, write_version_number


implicit none
private
      Integer,Parameter :: kind_type = selected_real_kind(15,307)

      Real, Allocatable, Dimension(:,:,:,:) :: TRNS

      Real,Private, Dimension(109)     :: PA
      Real,Private, Dimension(109,109) :: TRANSA
      Real(kind_type),Private, Dimension(109)   :: XA,CA,ETA,SEXPV
      Real(kind_type),Private                   :: CORE,UEXP,SEXP

      Real(kind_type),Private :: ZERO=0.0, ONE=1.0, TWO=2.0
!-----------------------------------------------------------------------

      Public   co2int, co2int_init, co2int_end, TRNS
      Private  RCTRNS,COEINT,QUADSR,SINTR2,Qintrp,PATH

!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: co2int.F90,v 13.0 2006/03/28 21:09:25 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------

      CONTAINS

!#######################################################################

      Subroutine co2int (nlev,ir,npurp,nkkk,unit1,unit2,ratio)

      Implicit None
!-----------------------------------------------------------------------
!
!      ------------   FUNCTION INTERPOLATER ROUTINE  ------------
!
!-----------------------------------------------------------------------
      Integer, Intent(IN) :: nlev,ir,npurp,nkkk,unit1,unit2
      Real,    Intent(IN) :: ratio
!-----------------------------------------------------------------------

      Real, Dimension(nlev+1) :: p,pd

      Real     P1,P2,TRNSLO,FACT15,FACT30,ratsm,ratstd
      Integer  i,j,kk,kkk,N,NLP1,NLP2,nmeth,ncalcs
!-----------------------------------------------------------------------

      NLP1=nlev+1
      NLP2=nlev+2

!-----------------------------------------------------------------------

!------ THE FOLLOWING ARE THE INPUT FORMATS -----

 100  FORMAT (4F20.14)

!-----------------------------------------------------------------------

!     CALCULATION OF PA -THE "TABLE" OF 109 GRID PRESSURES
!     NOTE-THIS CODE MUST NOT BE CHANGED BY THE USER!!!!!!!!!

      PA(1)=0.
      FACT15=10.**(1./15.)
      FACT30=10.**(1./30.)
      PA(2)=1.0E-3
      Do i=2,76
         PA(i+1)=PA(i)*FACT15
      EndDo
      Do i=77,108
         PA(i+1)=PA(i)*FACT30
      EndDo

      If (npurp == 1) Then
         ncalcs=4
      Else If (npurp == 2) Then
         ncalcs=2
      Else If (npurp == 3) Then
         ncalcs=2
      EndIf

!-------- allocate output transmission function arrays --------
      If (Allocated(TRNS)) DeAllocate (TRNS)
                             Allocate (TRNS(NLP1,NLP1,ncalcs,nkkk))

!=======================================================================
!***do loop on no. of temp profiles for each output calc (controlled
!   by nkkk,a function of freq. range)

      Do 410 kkk=1,nkkk
!-----------------------------------------------------------------------

!***read input lbl transmission fctn tapes (the no. depends on the
!   co2 amount and is controlled by unit1, unit2, a function of ratio)

      If (unit2 == 0) Then
         read (unit1,100) ((transa(i,j),i=1,109),j=1,109)
      Else
         If (ratio > 0.5 .and. ratio < 1.0) Then
             ratsm=0.5
             ratstd=1.0
         EndIf
         If (ratio > 1.0 .and. ratio < 2.0) Then
             ratsm=1.0
             ratstd=2.0
         EndIf
         If (ratio > 2.0 .and. ratio < 4.0) Then
             ratsm=2.0
             ratstd=4.0
         EndIf
         CALL RCTRNS (unit1,unit2,RATSTD,RATSM,RATIO,IR)
      EndIf

!***define interpolation coefficients in coeint

      Do i=1,109
         TRANSA(i,i)=1.0
      EndDo
      CALL COEINT (RATIO,IR)

!=======================================================================
!***do loop on number of output calculations (controlled by ncalcs,
!   a function of calculation purpose)

      Do 400 kk=1,ncalcs
!-----------------------------------------------------------------------

!---initialize transmission fctn array  

      Do i=1,NLP1
      Do j=1,NLP1
         TRNS(j,i,kk,kkk)=1.00
      EndDo
      EndDo

!***define pressure arrays pd,p according to the calc. no. kk and
!   npurp.define nmeth (interp method to be used)

      If (kk == 1) Then
         If (npurp == 1 .or. npurp == 3) Then
            pd=pd1013; p=plm1013
            nmeth=2 
         Else
            p=plm1013
            nmeth=1
         EndIf
      EndIf 

      If (kk == 2) Then
         If (npurp == 1) Then
            p=plm1013
            nmeth=1
         EndIf
         If (npurp == 2) Then
            p=plm810
            nmeth=1
         EndIf
         If (npurp == 3) Then
            pd=pd810; p=plm810
            nmeth=2
         EndIf
      EndIf

      If (kk == 3) Then
         pd=pd810; p=plm810
         nmeth=2
      EndIf

      If (kk == 4) Then
         p=plm810
         nmeth=1
      EndIf

      If (nmeth == 1) Then
         Do i=1,NLP1
         Do j=1,i
            IF (i == j) CYCLE
            P1=P(j)
            P2=P(i)
            CALL SINTR2 (P1,P2,TRNSLO)
            TRNS(j,i,kk,kkk)=TRNSLO
         EndDo
         EndDo
         Do i=1,NLP1
         Do j=i,NLP1
            TRNS(j,i,kk,kkk)=TRNS(i,j,kk,kkk)
         EndDo
         EndDo
      Else

!   perform method 1(point) calculations for (1,i) array elements. this
!   element will be used for cts calculations and is unneeded for other
!   aspects of radiation calcs, since we assume isothermal conditions
!   above the top pressure level.

         TRNS(1,1,kk,kkk)=1.00
         Do j=2,NLP1
            P1=0.
            P2=P(j)
            CALL SINTR2 (P1,P2,TRNSLO)
            TRNS(1,j,kk,kkk)=TRNSLO
         EndDo
         Do j=1,NLP1
         Do i=2,NLP1
            N=25
            IF (i /= j) N=3
            call quadsr (N,p(j),pd(i),pd(i-1),TRNS(i,j,kk,kkk))
         EndDo
         EndDo
      EndIf

!-----------------------------------------------------------------------
 400  Continue
 410  Continue
!=======================================================================

      End Subroutine co2int

!#######################################################################

      Subroutine co2int_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      End Subroutine co2int_init

!#######################################################################

      Subroutine co2int_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      End Subroutine co2int_end

!#######################################################################

      SUBROUTINE COEINT (RAT,IR)

      Implicit None
!-----------------------------------------------------------------------
!
!
!      THE TRANSMISSION FUNCTION BETWEEN P1 AND P2 IS ASSUMED TO  HAVE
!  THE FUNCTIONAL FORM
!               TAU(P1,P2)= 1.0-SQRT(C*LOG(1.0+X*PATH)),
!         WHERE
!               PATH(P1,P2)=((P1-P2)**2)*(P1+P2+CORE)/
!                           (ETA*(P1+P2+CORE)+(P1-P2))
!
!
!  THE PARAMETERS C AND X ARE FUNCTIONS OF P2, AND ARE TO BE DETERMINED,
!  WHILE CORE IS A PRESPECIFIED NUMBER.ETA IS A FUNCTION OF THE THEIR
!  PRODUCT (CX);IT IS OBTAITED ITERATIVELY. THE DERIVATION OF ALL THESE
!  VALUES WILL BE EXPLAINED IN A FORTHCOMING PAPER.
!      SUBROUTINE COEINT DETERMINES C(I) AND X(I) BY USING THE ACTUAL
!  VALUES OF TAU(P(I-2),P(I)) AND TAU(P(I-1),P(I)) AND THE PREVIOUS
!  ITERATION VALUE OF ETA.
!       DEFINE: 
!         PATHA=PATH(P(I),P(I-2),CORE,ETA)
!          PATHB=PATH(P(I),P(I-1),CORE,ETA);
!  THEN
!          R=(1-TAU(P(I),P(I-2)))/(1-TAU(P(I),P(I-1)))
!           = SQRT(LOG(1+X*PATHA)/LOG(1+X*PATHB)),
!  SO THAT
!          R**2= LOG(1+X*PATHA)/LOG(1+X*PATHB).
!  THIS EQUATION CAN BE SOLVED BY NEWTON'S METHOD FOR X AND THEN THE
!  RESULT USED TO FIND C. THIS IS REPEATED FOR EACH VALUE OF I GREATER
!  THAN 2 TO GIVE THE ARRAYS X(I) AND C(I).
!       NEWTON'S METHOD FOR SOLVING THE EQUATION
!           F(X)=0
!  MAKES USE OF THE LOOP XNEW= XOLD-F(XOLD)/F'(XOLD).
!  THIS IS ITERATED 20 TIMES, WHICH IS PROBABLY EXCESSIVE.
!  THE FIRST GUESS FOR ETA IS 3.2E-4*EXP(-P(I)/1000),WHICH HAS
!  BEEN FOUND TO BE FAIRLY REALISTIC BY EXPERIMENT; WE ITERATE 5 TIMES
!  (AGAIN,PROBABLY EXCESSIVELY) TO OBTAIN THE VALUES FOR C,X,ETA TO BE
!  USED FOR INTERPOLATION.
!     THERE ARE SEVERAL POSSIBLE PITFALLS: 
!        1) IN THE COURSE OF ITERATION, X MAY REACH A VALUE WHICH MAKES
!           1+X*PATHA NEGATIVE; IN THIS CASE THE ITERATION IS STOPPED,
!           AND AN ERROR MESSAGE IS PRINTED OUT.
!        2) EVEN IF (1) DOES NOT OCCUR, IT IS STILL POSSIBLE THAT X MAY
!           BE NEGATIVE AND LARGE ENOUGH TO MAKE
!           1+X*PATH(P(I),0,CORE,ETA) NEGATIVE. THIS IS CHECKED
!           FOR IN A FINAL LOOP, AND IF TRUE, A WARNING IS PRINTED OUT.
!
!-----------------------------------------------------------------------
      Real,    Intent(IN)  :: RAT
      Integer, Intent(IN)  :: IR
!-----------------------------------------------------------------------

      Real  PA2

      Real(kind_type) :: padi,padim,padim2,PATHA,PATHB,P0,R,REXP,XX,  &
                      ftest1,ftest2,xxlog,F1,F2,F,FPRIME,CHECK,drat
      Real(kind_type) :: PATH0(109),ETAP(109)
      Real(kind_type) :: SINV(4)
      Real(kind_type) :: small

      Integer  i,NP,LL

      Character(len=40) :: err_string

      DATA SINV /2.74992,2.12731,4.38111,.0832926/

!-----------------------------------------------------------------------

      CORE=5.000
      UEXP=0.90
      P0=0.7
      small = epsilon(CORE)

      Do i=1,109
         PA2=PA(i)*PA(i)
         SEXPV(i)=.505+2.0e-5*PA(i)+.035*(PA2-.25)/(PA2+.25)
      EndDo

      Do i=1,109
         ETA(i)=3.2e-4*EXP(-PA(i)/500.)
         ETAP(i)=ETA(i)
      EndDo

      Do NP=1,10

         Do i=3,109
            padi=pa(i)
            padim=pa(i-1)
            padim2=pa(i-2)
            SEXP=SEXPV(i)
            R=(one-TRANSA(i,i-2))/(one-TRANSA(i,i-1))
            REXP=R**(UEXP/SEXP)
            PATHA=(PATH(padi,padim2,CORE,ETA(i)))**UEXP
            PATHB=(PATH(padi,padim,CORE,ETA(i)))**UEXP
            XX=two*(PATHB*REXP-PATHA)/(PATHB*PATHB*REXP-PATHA*PATHA)
            Do LL=1,20
               ftest1=xx*patha
               ftest2=xx*pathb
!*** end iteration and solve if ftest1 is small or ftest2 is large
               If (ftest1 <= small) Then
                  xx=one
                  xa(i)=xx
                  ca(i)=(ONE-transa(i,i-2))**(uexp/sexp)/patha
                  GoTo 1011
               EndIf
               If (ftest2 >= 1.0e8) Then
                  xxlog=(Log(patha)-rexp*Log(pathb))/(rexp-ONE)
                  xx=Log(xxlog)
                  xa(i)=xx
                  ca(i)=(ONE-transa(i,i-2))**(uexp/sexp)/(xxlog+Log(patha))
                  GoTo 1011
               EndIf
               F1=Log(ONE+XX*PATHA)
               F2=Log(ONE+XX*PATHB)
               F=F1/F2-REXP
               FPRIME=(F2*PATHA/(ONE+XX*PATHA)-F1*PATHB/(ONE+XX*PATHB))/  &
                          (F2*F2)
               XX=XX-F/FPRIME
               CHECK=ONE+XX*PATHA
               If (CHECK <= 0.0) Then
                  Write  (err_string(1:37),360) i,LL,CHECK
  360             Format ('i=',i3,'LL=',i3,'CHECK=',f20.10)
                  CALL ERROR_MESG ('COEINT in CO2INT_MOD', &
                                    err_string(1:37), FATAL)
               EndIf
            EndDo

            CA(i)=(ONE-TRANSA(i,i-2))**(UEXP/SEXP)/  &
                    (Log(ONE+XX*PATHA)+small)
            XA(i)=XX
 1011    continue
         EndDo

         XA(2)=XA(3)
         XA(1)=XA(3)
         CA(2)=CA(3)
         CA(1)=CA(3)

         Do i=3,109
            padi=pa(i)
            PATH0(i)=(PATH(padi,ZERO,CORE,ETA(i)))**UEXP
            PATH0(i)=ONE+XA(i)*PATH0(i)
            If (PATH0(i) < 0.) then
               Write  (err_string(1:37),361) i
               CALL ERROR_MESG ('COEINT in CO2INT_MOD',   &
                                 err_string(1:37), WARNING)
            Endif
!del        If (PATH0(i) < 0.) Write (*,361) i,PATH0(i),XA(i)
         EndDo

         Do i=1,109
            drat=rat
            SEXP=SEXPV(i)
            ETAP(i)=ETA(i)
            ETA(i)=(SINV(IR)/drat)**(1./SEXP)*(CA(i)*XA(i))**(1./UEXP)
         EndDo

!-----------------------------------------------------------------------
!     THE ETA FORMULATION IS DETAILED IN SCHWARZKOPF AND FELS(1985).
!        THE QUANTITY SINV=(G*DELTANU)/(RCO2*D*S)
!      IN CGS UNITS,WITH D,THE DIFFUSICITY FACTOR=2, AND
!      S,THE SUM OF CO2 LINE STRENGTHS OVER THE 15UM CO2 BAND
!       ALSO,THE DENOMINATOR IS MULTIPLIED BY
!      1000 TO PERMIT USE OF MB UNITS FOR PRESSURE.
!        S IS ACTUALLY WEIGHTED BY B(250) AT 10 CM-1 WIDE INTERVALS,IN
!      ORDER TO BE CONSISTENT WITH THE METHODS USED TO OBTAIN THE LBL
!      1-BAND CONSOLIDATED TRANCMISSION FUNCTIONS.
!      FOR THE 490-850 INTERVAL (DELTANU=360,IR=1) SINV=2.74992.
!      (SLIGHTLY DIFFERENT FROM 2.7528 USED IN EARLIER VERSIONS)
!      FOR THE 490-670 INTERVAL (IR=2) SINV=2.12731
!      FOR THE 670-850 INTERVAL (IR=3) SINV=4.38111
!      FOR THE 2270-2380 INTERVAL (IR=4) SINV=0.0832926
!      SINV HAS BEEN OBTAINED USING THE 1982 AFGL CATALOG FOR CO2
!        RAT IS THE ACTUAL CO2 MIXING RATIO IN UNITS OF 330 PPMV,
!      LETTING USE OF THIS FORMULATION FOR ANY CO2 CONCENTRATION.
!-----------------------------------------------------------------------

!        Write  (*,366) (NP,i,CA(i),XA(i),ETA(i),SEXPV(i),i=1,109)
!366     Format (2i4,4e20.12)

      EndDo

 361  Format ('1+XA*PATH(PA(i),0) IS NEGATIVE,i= ',i3)
!361  Format (' **WARNING:** 1+XA*PATH(PA(i),0) IS NEGATIVE,i= ',i3,  &
!               /,20X,'PATH0(i)=',f16.6,' XA(i)=',f16.6)
!-----------------------------------------------------------------------

      END SUBROUTINE COEINT

!#######################################################################

      SUBROUTINE RCTRNS (ITAP1,ITAP2,RATSTD,RATSM,RATIO,IR)

      Implicit None
!-----------------------------------------------------------------------
!      RATSTD=VALUE OF HIGHER STD CO2 CONCENTRATION
!      RATSM=VALUE OF LOWER STD CO2 CONCENTRATION
!      RATIO=ACTUAL CO2 CONCENTRATION
!      THE 3 ABOVE QUANTITIES ARE IN UNITS OF 330 PPMV.
!-----------------------------------------------------------------------
      Integer, Intent(IN)  :: ITAP1,ITAP2
      Real,    Intent(IN)  :: RATSTD,RATSM,RATIO
      Integer, Intent(IN)  :: IR
!-----------------------------------------------------------------------

      Real :: TRNS1(109,109),TRNS2(109,109)

      Real      P1,P2,TRNSLO,TRNSPR,TRNSPM
      Integer   i,j
!-----------------------------------------------------------------------

!   READ IN TFS OF LOWER STD CO2 CONCENTRATION

      READ   (ITAP1,100) ((TRNS1(i,j),i=1,109),j=1,109)
100   FORMAT (4F20.14)

!   READ IN TFS OF HIGHER STD CO2 CONCENTRATION

      READ (ITAP2,100) ((TRANSA(i,j),i=1,109),j=1,109)

!-----------------------------------------------------------------------

!     CALL COEINT (RATSTD)
      CALL COEINT (RATSTD,IR)

      Do i=1,109
      Do j=1,i
         If (j == i) CYCLE

!  USING HIGHER CO2 CONCENTRATION,COMPUTE 1ST GUESS CO2 TFS FOR
!  ACTUAL CO2 CONCENTRATION.

         P2=(RATIO+RATSTD)*PA(i)/(2.*RATSTD) +  &
            (RATSTD-RATIO)*PA(j)/(2.*RATSTD)
         P1=(RATSTD-RATIO)*PA(i)/(2.*RATSTD) +  &
            (RATIO+RATSTD)*PA(j)/(2.*RATSTD)
         CALL SINTR2 (P1,P2,TRNSLO)
         TRNSPR=TRNSLO

!  USING HIGHER CO2 CONCENTRATION,COMPUTE 1ST GUESS CO2 TFS FOR
!  LOWER STD CO2 CONCENTRATION

         P2=(RATSM+RATSTD)*PA(i)/(2.*RATSTD) +  &
            (RATSTD-RATSM)*PA(j)/(2.*RATSTD)
         P1=(RATSTD-RATSM)*PA(i)/(2.*RATSTD) +  &
            (RATSM+RATSTD)*PA(j)/(2.*RATSTD)
         CALL SINTR2 (P1,P2,TRNSLO)
         TRNSPM=TRNSLO

!  COMPUTE TFS FOR CO2 CONCENTRATION GIVEN BY (RATIO).
!   STORE TEMPORARILY IN (TRNS2)

         TRNS2(j,i)=TRNSPR+(RATSTD-RATIO)*(TRNS1(j,i)-  &
          TRNSPM)/(RATSTD-RATSM)
         TRNS2(i,j)=TRNS2(j,i)

! WE NOW CAN OVERWRITE (TRNS1) AND STORE IN (TRNS1) THE 1ST GUESS
!  CO2 TFS FOR LOWER STD CO2 CONCENTRATION

         TRNS1(j,i)=TRNSLO
         TRNS1(i,j)=TRNSLO

      EndDo
      EndDo

!  SET DIAGONAL VALUES OF CO2 TFS TO UNITY

      Do i=1,109
         TRNS1(i,i)=1.0
         TRNS2(i,i)=1.0
      EndDo

!  NOW OUTPUT THE COMPUTED CO2 TFS FOR (RATIO) CO2 CONC. IN (TRANSA)

      Do i=1,109
      Do j=1,109
         TRANSA(j,i)=TRNS2(j,i)
      EndDo
      EndDo

!-----------------------------------------------------------------------

      END SUBROUTINE RCTRNS

!#######################################################################

      SUBROUTINE QUADSR (N,P,PD1,PD2,TRNS)

      Implicit None
!-----------------------------------------------------------------------
      Integer, Intent(IN)  :: N
      Real,    Intent(IN)  :: P,PD1,PD2
      Real,    Intent(OUT) :: TRNS
!-----------------------------------------------------------------------
!  Note:  PD1=PD(IA), PD2=PD(IA-1), P=P(JA)

      Real     WT(101)
      Real     TRNSNB,DP,PFIX,PVARY,P1,P2,TRNSLO
      Integer  i,kk,N2,N2P

      N2=2*N
      N2P=2*N+1

!------- WEIGHTS ARE CALCULATED ------
      WT(1)=1.
      Do i=1,N
         WT(2*i)=4.
         WT(2*i+1)=1.
      EndDo

      If (N > 1) Then
         Do i=2,N
            WT(2*i-1)=2.
         EndDo
      EndIf

      TRNSNB=0.
!!!!  DP=(PD(IA)-PD(IA-1))/N2
      DP=(PD1-PD2)/N2
      PFIX=P

      Do kk=1,N2P
!!!!     PVARY=PD(IA-1)+(kk-1)*DP
         PVARY=PD2+(kk-1)*DP
         IF (PVARY >= PFIX) P2=PVARY
         IF (PVARY >= PFIX) P1=PFIX
         IF (PVARY <  PFIX) P1=PVARY
         IF (PVARY <  PFIX) P2=PFIX
         CALL SINTR2 (P1,P2,TRNSLO)
         TRNSNB=TRNSNB+TRNSLO*WT(kk)
      EndDo

!!!!  TRNS(IA,JA)=TRNSNB*DP/(3.*(PD(IA)-PD(IA-1)))
      TRNS       =TRNSNB*DP/(3.*(PD1-PD2))

!-----------------------------------------------------------------------

      END SUBROUTINE QUADSR

!#######################################################################

      SUBROUTINE SINTR2 (P1,P2,TRNSLO)

      Implicit None
!-----------------------------------------------------------------------
      Real,Intent(IN)  :: P1,P2
      Real,Intent(OUT) :: TRNSLO
!-----------------------------------------------------------------------

      Real(kind_type) :: &
      p1d,p2d,padi,padip,padj,padjp,padip2,padjp2,PETA,paieta,         &
      paiet1,ETAP,PIPMPI,UP2P1,suexp,TRIP,TRI,TIJ,TIPJ,TIJP,TIPJP,     &
      UIJ,UIPJ,UIJP,UIPJP,PRODI,PRODIP,PROD,XINT,CINT,AIJ,AIJP,AIPJ,   &
      AIPJP,EIJ,EIPJ,EIJP,EIPJP,DTDJ,DTDPJ,EPIP1,EPIPP1,EPP2P1,TIP2J,  &
      TIP2JP,TI2J2,TIJP2,TIPJP2,UIP2J,UIJP2,UIPJP2,UI2J2,UIP2JP,       &
      AIJP2,AIPJP2,AIP2J,AIP2JP,AI2J2,EIP2J,EIP2JP,EIJP2,EIPJP2,       &
      EI2J2,EI,EP,EP2,EPSIL

      Integer   i,j,k,ieta

!---find indices for pa corresponding to p1 and p2 (to use for 
!   pressure interpolation

      If (p2 <= pa(1)) Then
         i=1
      EndIf

      Do k=1,108
        If (p2 > pa(k).and.p2 <= pa(k+1)) Then
           i=k
        EndIf
      EndDo

      If (p2 > pa(109)) Then
        i=108
      EndIf

      If (p1 <= pa(1)) Then
         j=1
      EndIf

      Do k=1,108
        If (p1 > pa(k) .and. p1 <= pa(k+1)) Then
           j=k
        EndIf
      EndDo

      If (p1 > pa(109)) Then
        j=108
      EndIf

!--define real(kind_type) quantities for pressures used in calc.

      p1d=p1
      p2d=p2
      padi=pa(i)
      padip=pa(i+1)
      padj=pa(j)
      padjp=pa(j+1)
      If (i < 108) padip2=pa(i+2)
      If (j < 108) padjp2=pa(j+2)

!  DETERMINE ETAP,THE VALUE OF ETA TO USE BY LINEAR INTERPOLATION
!    FOR PETA(=0.5*(P1+P2))
      PETA=p2d 

!---if peta=p2d,ieta will equal i
      ieta=i
      paieta=pa(ieta)
      paiet1=pa(ieta+1)
      ETAP=ETA(IETA)+(p2d-paieta)*(ETA(ieta+1)-ETA(IETA))/  &
             (paiet1-paieta)
      SEXP=SEXPV(IETA)+(p2d-paieta)*(SEXPV(ieta+1)-SEXPV(IETA))/  &
             (paiet1-paieta)
      PIPMPI=padip-padi
      UP2P1=(PATH(p2d,p1d,CORE,ETAP))**UEXP
      suexp=sexp/uexp

      If (i <= j) Then
        TRIP=(CA(i+1)*Log(ONE+XA(i+1)*UP2P1))**suexp
        TRI=(CA(i)*Log(ONE+XA(i)*UP2P1))**suexp
        TRNSLO=ONE-((padip-p2d)*TRI+(p2d-padi)*TRIP)/PIPMPI
      EndIf

      If (i > j) Then
         TIJ=TRANSA(i,j)
         TIPJ=TRANSA(i+1,j)
         TIJP=TRANSA(i,j+1)
         TIPJP=TRANSA(i+1,j+1)
         UIJ=(PATH(padi,padj,CORE,ETAP))**UEXP
         UIPJ=(PATH(padip,padj,CORE,ETAP))**UEXP
         UIJP=(PATH(padi,padjp,CORE,ETAP))**UEXP
         UIPJP=(PATH(padip,padjp,CORE,ETAP))**UEXP
         PRODI=CA(i)*XA(i)
         PRODIP=CA(i+1)*XA(i+1)
         PROD=((padip-p2d)*PRODI+(p2d-padi)*PRODIP)/PIPMPI
         XINT=((padip-p2d)*XA(i)+(p2d-padi)*XA(i+1))/PIPMPI
         CINT=PROD/XINT
         AIJ=(CINT*Log(ONE+XINT*UIJ))**suexp
         AIJP=(CINT*Log(ONE+XINT*UIJP))**suexp
         AIPJ=(CINT*Log(ONE+XINT*UIPJ))**suexp
         AIPJP=(CINT*Log(ONE+XINT*UIPJP))**suexp
         EIJ=TIJ+AIJ
         EIPJ=TIPJ+AIPJ
         EIJP=TIJP+AIJP
         EIPJP=TIPJP+AIPJP
         DTDJ=(EIJP-EIJ)/(padjp-padj)
         DTDPJ=(EIPJP-EIPJ)/(padjp-padj)
         EPIP1=EIJ+DTDJ*(p1d-padj)
         EPIPP1=EIPJ+DTDPJ*(p1d-padj)
         EPP2P1=((padip-p2d)*EPIP1+(p2d-padi)*EPIPP1)/PIPMPI
         TRNSLO=EPP2P1-(CINT*Log(ONE+XINT*UP2P1))**suexp
      EndIf

      If (i < 108 .and. j < 108 .and. i > j+2) Then
         TIP2J=TRANSA(i+2,j)
         TIP2JP=TRANSA(i+2,j+1)
         TI2J2=TRANSA(i+2,j+2)
         TIJP2=TRANSA(i,j+2)
         TIPJP2=TRANSA(i+1,j+2)
         UIP2J=(PATH(padip2,padj,CORE,ETAP))**UEXP
         UIJP2=(PATH(padi,padjp2,CORE,ETAP))**UEXP
         UIPJP2=(PATH(padip,padjp2,CORE,ETAP))**UEXP
         UI2J2=(PATH(padip2,padjp2,CORE,ETAP))**UEXP
         UIP2JP=(PATH(padip2,padjp,CORE,ETAP))**UEXP
         AIJP2=(CINT*Log(ONE+XINT*UIJP2))**suexp
         AIPJP2=(CINT*Log(ONE+XINT*UIPJP2))**suexp
         AIP2J=(CINT*Log(ONE+XINT*UIP2J))**suexp
         AIP2JP=(CINT*Log(ONE+XINT*UIP2JP))**suexp
         AI2J2=(CINT*Log(ONE+XINT*UI2J2))**suexp
         EIP2J=TIP2J+AIP2J
         EIP2JP=TIP2JP+AIP2JP
         EIJP2=TIJP2+AIJP2
         EIPJP2=TIPJP2+AIPJP2
         EI2J2=TI2J2+AI2J2
         CALL QINTRP(padj,padjp,padjp2,EIJ,EIJP,EIJP2,p1d,EI)
         CALL QINTRP(padj,padjp,padjp2,EIPJ,EIPJP,EIPJP2,p1d,EP)
         CALL QINTRP(padj,padjp,padjp2,EIP2J,EIP2JP,EI2J2,p1d,EP2)
         CALL QINTRP(padi,padip,padip2,EI,EP,EP2,p2d,EPSIL)
         TRNSLO=EPSIL-(CINT*Log(ONE+XINT*UP2P1))**suexp
      EndIf

!-----------------------------------------------------------------------

      END SUBROUTINE SINTR2

!#######################################################################

      Subroutine Qintrp (XM,X0,XP,FM,F0,FP,X,F)

      Implicit None
!-----------------------------------------------------------------------
      Real(kind_type), Intent(IN)  :: XM,X0,XP,FM,F0,FP,X
      Real(kind_type), Intent(OUT) :: F
!-----------------------------------------------------------------------
      Real(kind_type) :: D1,D2,B,A,DEL

      D1=(FP-F0)/(XP-X0)
      D2=(FM-F0)/(XM-X0)
      B=(D1-D2)/(XP-XM)
      A=D1-B*(XP-X0)
      DEL=(X-X0)
      F=F0+DEL*(A+DEL*B)

!-----------------------------------------------------------------------

      End Subroutine Qintrp

!#######################################################################

      FUNCTION PATH (A,B,C,E)

      Implicit None
      Real(kind_type), Intent(IN) :: A,B,C,E

      Real(kind_type) :: PEXP
      Real(kind_type) :: PATH

      PEXP=1./SEXP
      PATH=((A-B)**PEXP*(A+B+C))/(E*(A+B+C)+(A-B)**(PEXP-1.))

      END FUNCTION PATH

!#######################################################################

      End Module co2int_mod

