
      Module CO2_Data_Mod

!-----------------------------------------------------------------------

      use fs_profile_mod, ONLY:  fs_profile
      Use     co2int_mod, ONLY:  co2int, TRNS
      Use        fms_mod, ONLY:  open_namelist_file, mpp_pe,  &
                                 Error_Mesg, FATAL, close_file,  &
                                 write_version_number, mpp_root_pe, open_file

implicit none
private

!-----------------------------------------------------------------------
!
!   Pretabulated co2 transmission functions, evaluated using the
!   methods of Fels and Schwarzkopf (1981) and Schwarzkopf and
!   Fels (1985). 
!
!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 560-800 cm-1 band. also included are the
!   standard temperatures and the weighting function.
!   This data was formally in COMMON /CO2BD3/.
!
!       CO251    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO258    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) 
!                     WITH P(SFC)= ^810 MB
!       CDT51    =  FIRST TEMPERATURE DERIVATIVE OF CO251 
!       CDT58    =  FIRST TEMPERATURE DERIVATIVE OF CO258 
!       C2D51    =  SECOND TEMPERATURE DERIVATIVE OF CO251
!       C2D58    =  SECOND TEMPERATURE DERIVATIVE OF CO251
!       CO2M51   =  TRANSMISSION FCTNS FOR T0 FOR ADJACENT PRESSURE 
!                      LEVELS, WITH NO PRESSURE QUADRATURE. USED FOR
!                      NEARBY LAYER COMPUTATIONS. P(SFC)=1013.25 MB 
!       CO2M58   =  SAME AS CO2M51,WITH P(SFC)= ^810 MB 
!       CDTM51   =  FIRST TEMPERATURE DERIVATIVE OF CO2M51
!       CDTM58   =  FIRST TEMPERATURE DERIVATIVE OF CO2M58
!       C2DM51   =  SECOND TEMPERATURE DERIVATIVE OF CO2M51 
!       C2DM58   =  SECOND TEMPERATURE DERIVATIVE OF CO2M58 
!       STEMP    =  STANDARD TEMPERATURES FOR MODEL PRESSURE LEVEL
!                      STRUCTURE WITH P(SFC)=1013.25 MB 
!       GTEMP    =  WEIGHTING FUNCTION FOR MODEL PRESSURE LEVEL 
!                      STRUCTURE WITH P(SFC)=1013.25 MB.
!       B0       =  TEMP. COEFFICIENT USED FOR CO2 TRANS. FCTN. 
!                      CORRECTION FOR T(K). (SEE REF. 4 AND BD3)
!       B1       =  TEMP. COEFFICIENT, USED ALONG WITH B0 
!       B2       =  TEMP. COEFFICIENT, USED ALONG WITH B0 
!       B3       =  TEMP. COEFFICIENT, USED ALONG WITH B0 

      Real, Allocatable, Dimension(:,:) :: CO251,CO258,CDT51,CDT58
      Real, Allocatable, Dimension(:,:) :: C2D51,C2D58
      Real, Allocatable, Dimension(:)   :: CO2M51,CO2M58,CDTM51,CDTM58
      Real, Allocatable, Dimension(:)   :: C2DM51,C2DM58
      Real, Allocatable, Dimension(:)   :: STEMP,GTEMP
      Real                              :: B0,B1,B2,B3

!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 560-670 cm-1 part of the 15 um co2 band. 
!   This data was formally in COMMON /CO2BD2/.
!
!       CO231    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO238    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB
!       CDT31    =  FIRST TEMPERATURE DERIVATIVE OF CO231
!       CDT38    =  FIRST TEMPERATURE DERIVATIVE OF CO238
!       C2D31    =  SECOND TEMPERATURE DERIVATIVE OF CO231
!       C2D38    =  SECOND TEMPERATURE DERIVATIVE OF CO231

      Real, Allocatable, Dimension(:) :: CO231,CO238,CDT31,CDT38
      Real, Allocatable, Dimension(:) :: C2D31,C2D38

!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 670-800 part of the 15 um co2 band.
!   This data was formally in COMMON /CO2BD4/.
!
!       CO271    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO278    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB
!       CDT71    =  FIRST TEMPERATURE DERIVATIVE OF CO271
!       CDT78    =  FIRST TEMPERATURE DERIVATIVE OF CO278
!       C2D71    =  SECOND TEMPERATURE DERIVATIVE OF CO271
!       C2D78    =  SECOND TEMPERATURE DERIVATIVE OF CO271

      Real, Allocatable, Dimension(:) :: CO271,CO278,CDT71,CDT78
      Real, Allocatable, Dimension(:) :: C2D71,C2D78

!-----------------------------------------------------------------------
!
!   co2 transmission functions for the 2270-2380 part of the
!   4.3 um co2 band. THis data was formally in COMMON /CO2BD5/.
!
!       CO211    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO218    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB

      Real, Allocatable, Dimension(:) :: CO211,CO218

! 
!-----------------------------------------------------------------------

      Integer, Private :: LP1=0,LMAX=0

!-----------------------------------------------------------------------
!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: co2_data.F90,v 13.0 2006/03/28 21:09:22 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.

!-----------------------------------------------------------------------

      Public   CO2_Data, Write_CO2_Data, Read_CO2_Data, &
               co2_data_init, co2_data_end

      Public   CO251,  CO258,  CDT51,  CDT58,  C2D51,  C2D58,  &
               CO2M51, CO2M58, CDTM51, CDTM58, C2DM51, C2DM58, &
               STEMP,  GTEMP,  B0,     B1,     B2,     B3,     &
               CO231,  CO238,  CDT31,  CDT38,  C2D31,  C2D38,  &
               CO271,  CO278,  CDT71,  CDT78,  C2D71,  C2D78,  &
               CO211,  CO218

      CONTAINS

!#######################################################################

      Subroutine CO2_Data (co2std, ratio, Pref)

      Implicit None
!-----------------------------------------------------------------------
      Real, Intent(IN) :: co2std, ratio
      Real, Intent(IN) :: Pref(:,:)
!-----------------------------------------------------------------------
!   CO2STD = standard co2 vol. mixing ratio (either 300 or 330 ppmv)
!   RATIO  = co2 vol. mixing ratio in units of the standard vol. 
!            mixing ratio (must lie between 0.5 and 4.0)
!   PREF   = reference pressure levels
!-----------------------------------------------------------------------
      Integer    ir,iq,npurp,nkkk,unit1,unit2,m,n,ntap
      Real       co2mix,ccomp
      Real, Dimension(4) :: cstd = (/ 0.5, 1.0, 2.0, 4.0 /)
!-----------------------------------------------------------------------
      Character(len=14), Dimension(3) :: files =  &
         (/ 'INPUT/cns_300_', 'INPUT/cns_330_', 'INPUT/cns_600_' /)
      Character(len=6), Dimension(4) :: bands =  &
        (/ '490850', '490670', '670850', '43um  ' /)
!-----------------------------------------------------------------------


!----- check input values -----

      If (ratio < 0.5 .or. ratio > 4.0)  Call Error_Mesg  &
       ('CO2_Data', 'ratio > 4.0 or ratio < 0.5', FATAL)

      co2mix=co2std*ratio

!-----------------------------------------------------------------------

      B0 = -.51926410E-4
      B1 = -.18113332E-3
      B2 = -.10680132E-5
      B3 = -.67303519E-7

!-----------------------------------------------------------------------
!---------- has this data been previously allocated ? ------------------

      If (Size(Pref,1) .ne. LP1) Then
         LP1=Size(Pref,1); LMAX=LP1-1

         If (Allocated(CDT51)) Then
            DeAllocate (CDT51, CO251, C2D51, CDT58, CO258, C2D58)
            DeAllocate (CDT31, CO231, C2D31, CDT38, CO238, C2D38)
            DeAllocate (CDT71, CO271, C2D71, CDT78, CO278, C2D78)
            DeAllocate (CO211, CO218, STEMP, GTEMP)
            DeAllocate (CDTM51, CO2M51, C2DM51, CDTM58, CO2M58, C2DM58)
         EndIf

!  ----- For the 560-800 cm-1 bandwidth -----

         Allocate (CDT51(LP1,LP1), CO251(LP1,LP1), C2D51(LP1,LP1),  &
                   CDT58(LP1,LP1), CO258(LP1,LP1), C2D58(LP1,LP1))

!  ----- For the 560-670 cm-1 bandwidth -----

         Allocate (CDT31(LP1), CO231(LP1), C2D31(LP1),  &
                   CDT38(LP1), CO238(LP1), C2D38(LP1))

!  ----- For the 670-800 cm-1 bandwidth -----

         Allocate (CDT71(LP1), CO271(LP1), C2D71(LP1),  &
                   CDT78(LP1), CO278(LP1), C2D78(LP1))

!  ----- For the 2270-2380 cm-1 bandwidth -----

         Allocate (CO211(LP1), CO218(LP1))

! ----- For the 560-800 cm-1 bandwidth -----

         Allocate (CDTM51(LMAX), CO2M51(LMAX), C2DM51(LMAX),  &
                   CDTM58(LMAX), CO2M58(LMAX), C2DM58(LMAX))

!  STEMP IS THE US STANDARD ATMOSPHERES,1976,AT N18 PRESSURES
!  WHERE PSTAR=1013.25 MB
!  THE WEIGHTING FUNCTION GTEMP=P(K)**0.2*(1.+P(K)/30000.)**0.8/
!  1013250.,WHERE P(K)=PRESSURE,N18 DATA LEVELS FOR PSTAR=
!  1013250.

         Allocate (STEMP(LP1), GTEMP(LP1))

      EndIf

!-----------------------------------------------------------------------
!---- compute profiles -----

      Call fs_profile (Pref,STEMP,GTEMP)

!-----------------------------------------------------------------------
!
!   Use one tape (no interpolation) if the mixing ratio is
!   sufficiently close to the standard value, also do not interpolate
!   if the mixing ratio is sufficiently close to multiples
!   (0.5,2.0,4.0) of the standard value
!
!-----------------------------------------------------------------------

      If (co2std <= 310.) Then
          m=1
      Else if (co2std >  310. .and. co2std <= 350.) Then
          m=2
      Else
          m=3
      EndIf

      Do n=1,4
         ccomp=abs(co2mix-cstd(n)*co2std)
         if (ccomp <= 1.0e-4) then
            ntap=1
            EXIT
         else
            ntap=2
         endif
      EndDo

!-----------------------------------------------------------------------
!        ir=1:  lbl transmissions over 490-850 cm-1
!        ir=2:  lbl transmissions over 490-670 cm-1
!        ir=3:  lbl transmissions over 670-850 cm-1
!        ir=4:  lbl transmissions over 2270-2380 cm-1
!-----------------------------------------------------------------------

      Do ir = 1,4

!-----------------------------------------------------------------------
!
!   read in indices npurp and nkkk.
!   for GCM radiation codes, the values to be used for
!   npurp and nkkk vary with ir as follows:
!     ir=1:   npurp=1,nkkk=3
!     ir=2:   npurp=3,nkkk=3
!     ir=3:   npurp=3,nkkk=3
!     ir=4:   npurp=3,nkkk=1
!
!   read in npurp, an index giving the kinds of interps
!   desired, according to the following values:
!     npurp=1: purpose 1 and purpose 2 calcs. for p(sfc)
!              of 1013.25 and 810.6 mb.
!     npurp=2: purpose 1 calcs. only
!     npurp=3: purpose 2 calcs. only
!   purposes 1 and 2 are explained in comments for subroutine co2int
!
!   read in nkkk, an index giving the no. of temp.'
!   profiles to be calculated. use the following values:'
!     nkkk=1: calculate only the (T0) profile'
!     nkkk=3: calculate the (T0,T+ and T-) profiles'
!             (normal case).'
!      note: if ir=4, nkkk must be 1'
!
!-----------------------------------------------------------------------

         Select Case (ir)
            Case (1)
               npurp=1; nkkk=3
            Case (2:3)
               npurp=3; nkkk=3
            Case (4)
               npurp=3; nkkk=1
         End Select

!-----------------------------------------------------------------------
!
!***open lbl co2 transmission functions pior to executing the interpol-
!   ation pgm. this will be user-dependent. additionally, it depends
!   on ir (freq. range) and on ratio (co2 amt). the files below 
!   assume that ratio lies between 1 and 2, and that the file and
!   directory names are as specified below.
!
!-----------------------------------------------------------------------

         If (ntap == 1) Then
            unit1 = open_file (file=files(m)//bands(ir), action='read')
            unit2 = 0
         EndIf

         If (ntap == 2) Then
            unit1 = open_file (file=files(1)//bands(ir), action='read')
            unit2 = open_file (file=files(3)//bands(ir), action='read')
         EndIf

!-----------------------------------------------------------------------
!----- interpolate co2 transmission fctns -----

         Call co2int (LMAX,ir,npurp,nkkk,unit1,unit2,ratio)  

!RSH     call close_file (unit1, status='keep') 
!RSH     If (ntap == 2) call close_file (unit2, status='keep')
         call close_file (unit1               ) 
         If (ntap == 2) call close_file (unit2               )

!-----------------------------------------------------------------------
!     Load transmission functions into data arrays


         If (ir == 1) Then
            iq=ir
            Call co2ins (LMAX,1,3,iq)
            Call co2in1 (LMAX,2,4,iq)
         Else
            iq=ir
            If (ir == 4) iq=5
            Call co2ins(LMAX,1,2,iq)
         EndIf

!-----------------------------------------------------------------------
      EndDo
!-----------------------------------------------------------------------

      End Subroutine CO2_Data

!#######################################################################

      Subroutine co2ins (nlev,itin,itin1,iq)

      Implicit None
!-----------------------------------------------------------------------
      Integer, Intent(IN) :: nlev,itin,itin1,iq

      Real, Dimension(nlev+1,nlev+1) :: DCDT8,DCDT10,CO2PO,CO2800,  &
                                        CO2PO1,CO2801,CO2PO2,CO2802,  &
                                        D2CT8,D2CT10

      Integer   i,j,L,LP1,JMAX
      Real      C1,C2
!-----------------------------------------------------------------------

      L=nlev
      LP1=L+1

      CO2PO (:,:)=TRNS(:,:,itin ,1)
      CO2800(:,:)=TRNS(:,:,itin1,1)
      If (iq >= 1 .and. iq <= 4) Then
         CO2PO1(:,:)=TRNS(:,:,itin ,2)
         CO2801(:,:)=TRNS(:,:,itin1,2)
         CO2PO2(:,:)=TRNS(:,:,itin ,3)
         CO2802(:,:)=TRNS(:,:,itin1,3)
      EndIf

!-----------------------------------------------------------------------
!   THE FOLLOWING CODE IS REWRITTEN SO THAT THE RADIATIVE BANDS ARE: 
!
!        iq=1    560-800     (CONSOL.=490-850)
!        iq=2    560-670     (CONSOL.=490-670)
!        iq=3    670-800     (CONSOL.=670-850)
!        iq=4    560-760 (ORIGINAL CODE)   (CONSOL.=490-850)
!        iq=5   2270-2380     CONSOL=2270-2380
!
!  THE FOLLOWING LOOP OBTAINS TRANSMISSION FUNCTIONS FOR BANDS
!  USED IN RADIATIVE MODEL CALCULATIONS,WITH THE EQUIVALENT
!  WIDTHS KEPT FROM THE ORIGINAL CONSOLIDATED CO2 TF'S.
!      NOTE: ALTHOUGH THE BAND TRANSMISSION FUNCTIONS ARE
!  COMPUTED FOR ALL RADIATIVE BANDS, AS OF 9/28/88, THEY
!  ARE WRITTEN OUT IN FULL ONLY FOR THE FULL 15 UM BAND CASES
!  (iq=1,4).  IN OTHER CASES, THE TRANSMISSIVITIES (1,K) ARE
!  WRITTEN OUT, AS THESE ARE THE ONLY ONES NEEDED FOR CTS
!  CALCULATIONS.  ALSO, FOR THE 4.3 UM BAND (iq=5) THE TEMP.
!  DERIVATIVE TERMS ARE NOT WRITTEN OUT, AS THEY ARE UNUSED.
!-----------------------------------------------------------------------

      If (iq == 1) Then
         C1=1.5
         C2=0.5
         JMAX=LP1
      EndIf
      If (iq == 2) Then
        C1=18./11.
        C2=7./11.
        JMAX=1
      EndIf
      If (iq == 3) Then
        C1=18./13.
        C2=5./13.
        JMAX=1
      EndIf
      If (iq == 4) Then
        C1=1.8
        C2=0.8
        JMAX=LP1
      EndIf
      If (iq == 5) Then
        C1=1.0
        C2=0.0
        JMAX=1
      EndIf

      Do i=1,LP1
      Do j=1,LP1
         CO2PO(j,i)=C1*CO2PO(j,i)-C2
         CO2800(j,i)=C1*CO2800(j,i)-C2
      EndDo
      EndDo

      If (iq >= 1 .and. iq <= 4) Then
        Do i=1,LP1
        Do j=1,LP1
          CO2PO1(j,i)=C1*CO2PO1(j,i)-C2
          CO2801(j,i)=C1*CO2801(j,i)-C2
          CO2PO2(j,i)=C1*CO2PO2(j,i)-C2
          CO2802(j,i)=C1*CO2802(j,i)-C2
        EndDo
        EndDo
        Do j=1,LP1
        Do i=1,LP1
         DCDT8(i,j)=.02*(CO2801(i,j)-CO2802(i,j))*100.
         DCDT10(i,j)=.02*(CO2PO1(i,j)-CO2PO2(i,j))*100.
         D2CT8(i,j)=.0016*(CO2801(i,j)+CO2802(i,j)-2.*CO2800(i,j))*1000.
         D2CT10(i,j)=.0016*(CO2PO1(i,j)+CO2PO2(i,j)-2.*CO2PO(i,j))*1000.
        EndDo
        EndDo
      EndIf

!-----------------------------------------------------------------------
      If (iq == 1) Then
         CDT51=DCDT10
         CO251=CO2PO
         C2D51=D2CT10
         CDT58=DCDT8
         CO258=CO2800
         C2D58=D2CT8
      EndIf
!-----------------------------------------------------------------------
      If (iq == 2) Then
         CDT31=DCDT10(1,:)
         CO231=CO2PO (1,:)
         C2D31=D2CT10(1,:)
         CDT38=DCDT8 (1,:)
         CO238=CO2800(1,:)
         C2D38=D2CT8 (1,:)
      EndIf
!-----------------------------------------------------------------------
      If (iq == 3) Then
         CDT71=DCDT10(1,:)
         CO271=CO2PO (1,:)
         C2D71=D2CT10(1,:)
         CDT78=DCDT8 (1,:)
         CO278=CO2800(1,:)
         C2D78=D2CT8 (1,:)
      EndIf
!-----------------------------------------------------------------------
      If (iq == 4) Then
         Call Error_Mesg ('co2ins', 'iq cannot equal 4', FATAL)
      EndIf
!-----------------------------------------------------------------------
      If (iq == 5) Then
         CO211=CO2PO (1,:)
         CO218=CO2800(1,:)
      EndIf
!-----------------------------------------------------------------------

      End Subroutine co2ins

!#######################################################################

      Subroutine co2in1 (nlev,itin,itin1,iq)

      Implicit None
!-----------------------------------------------------------------------
!
!                       CO2INS FOR METHOD 1
!
!-----------------------------------------------------------------------

      Integer, Intent(IN)    :: nlev,itin,itin1,iq

      Real,   Dimension(nlev+1,nlev+1) :: DCDT8,DCDT10,CO2PO,CO2800,   &
                                          CO2PO1,CO2801,CO2PO2,CO2802, &
                                          D2CT8,D2CT10

      Integer  L,LP1,i,j
      Real     C1,C2
!-----------------------------------------------------------------------

      L=nlev
      LP1=L+1

      CO2PO (:,:)=TRNS(:,:,itin ,1)
      CO2800(:,:)=TRNS(:,:,itin1,1)
      If (IQ >= 1 .and. IQ <= 4) Then
         CO2PO1(:,:)=TRNS(:,:,itin ,2)
         CO2801(:,:)=TRNS(:,:,itin1,2)
         CO2PO2(:,:)=TRNS(:,:,itin ,3)
         CO2802(:,:)=TRNS(:,:,itin1,3)
      EndIf

!-----------------------------------------------------------------------
!   THE FOLLOWING CODE IS REWRITTEN SO THAT THE RADIATIVE BANDS ARE:
!
!        iq=1    560-800     (CONSOL.=490-850)
!        iq=2    560-670     (CONSOL.=490-670)
!        iq=3    670-800     (CONSOL.=670-850)
!        iq=4    560-760 (ORIGINAL CODE)   (CONSOL.=490-850)
!        iq=5   2270-2380     CONSOL=2270-2380
!
!  THE FOLLOWING LOOP OBTAINS TRANSMISSION FUNCTIONS FOR BANDS
!  USED IN RADIATIVE MODEL CALCULATIONS,WITH THE EQUIVALENT
!  WIDTHS KEPT FROM THE ORIGRNAL CONSOLIDATED CO2 TF'S.
!-----------------------------------------------------------------------

      If (iq == 1) Then
         C1=1.5
         C2=0.5
      EndIf
      If (iq == 2) Then
        C1=18./11.
        C2=7./11.
      EndIf
      If (iq == 3) Then
        C1=18./13.
        C2=5./13.
      EndIf
      If (iq == 4) Then
        C1=1.8
        C2=0.8
      EndIf
      If (iq == 5) Then
        C1=1.0
        C2=0.0
      EndIf

      Do i=1,LP1
      Do j=1,LP1
         CO2PO(j,i)=C1*CO2PO(j,i)-C2
         CO2800(j,i)=C1*CO2800(j,i)-C2
      EndDo
      EndDo

      If (iq >= 1 .and. iq <= 4) Then
        Do i=1,LP1
        Do j=1,LP1
         CO2PO1(j,i)=C1*CO2PO1(j,i)-C2
         CO2801(j,i)=C1*CO2801(j,i)-C2
         CO2PO2(j,i)=C1*CO2PO2(j,i)-C2
         CO2802(j,i)=C1*CO2802(j,i)-C2
        EndDo
        EndDo

        Do j=1,LP1
        Do i=1,LP1
         DCDT8(i,j)=.02*(CO2801(i,j)-CO2802(i,j))*100.
         DCDT10(i,j)=.02*(CO2PO1(i,j)-CO2PO2(i,j))*100.
         D2CT8(i,j)=.0016*(CO2801(i,j)+CO2802(i,j)-2.*CO2800(i,j))*1000.
         D2CT10(i,j)=.0016*(CO2PO1(i,j)+CO2PO2(i,j)-2.*CO2PO(i,j))*1000.
        EndDo
        EndDo
      EndIf

!-----------------------------------------------------------------------
      If (iq == 1) Then
         Do i=1,nlev
            CDTM51(i)=DCDT10(i,i+1)
            CO2M51(i)=CO2PO (i,i+1)
            C2DM51(i)=D2CT10(i,i+1)
            CDTM58(i)=DCDT8 (i,i+1)
            CO2M58(i)=CO2800(i,i+1)
            C2DM58(i)=D2CT8 (i,i+1)
         EndDo
      EndIf
!-----------------------------------------------------------------------

      End Subroutine co2in1

!#######################################################################

      Subroutine Write_CO2_Data

      Implicit None
!-----------------------------------------------------------------------
      Integer  i,k,nlev,nlevp1,unit
!-----------------------------------------------------------------------

      nlevp1=Size(CDT51,1); nlev=nlevp1-1

      unit = open_file (file='co2data', action='write')

      Do k=1,nlevp1; Write (unit,101) (CDT51(i,k),i=1,nlevp1); EndDo
      Do k=1,nlevp1; Write (unit,101) (CO251(i,k),i=1,nlevp1); EndDo
      Do k=1,nlevp1; Write (unit,101) (C2D51(i,k),i=1,nlevp1); EndDo
      Do k=1,nlevp1; Write (unit,101) (CDT58(i,k),i=1,nlevp1); EndDo
      Do k=1,nlevp1; Write (unit,101) (CO258(i,k),i=1,nlevp1); EndDo
      Do k=1,nlevp1; Write (unit,101) (C2D58(i,k),i=1,nlevp1); EndDo

                     Write (unit,101) (CDT31(i),i=1,nlevp1)
                     Write (unit,101) (CO231(i),i=1,nlevp1)
                     Write (unit,101) (C2D31(i),i=1,nlevp1)
                     Write (unit,101) (CDT38(i),i=1,nlevp1)
                     Write (unit,101) (CO238(i),i=1,nlevp1)
                     Write (unit,101) (C2D38(i),i=1,nlevp1)

                     Write (unit,101) (CDT71(i),i=1,nlevp1)
                     Write (unit,101) (CO271(i),i=1,nlevp1)
                     Write (unit,101) (C2D71(i),i=1,nlevp1)
                     Write (unit,101) (CDT78(i),i=1,nlevp1)
                     Write (unit,101) (CO278(i),i=1,nlevp1)
                     Write (unit,101) (C2D78(i),i=1,nlevp1)

                     Write (unit,101) (CO211(i),i=1,nlevp1)
                     Write (unit,101) (CO218(i),i=1,nlevp1)

                     Write (unit,101) (CDTM51(i),i=1,nlev)
                     Write (unit,101) (CO2M51(i),i=1,nlev)
                     Write (unit,101) (C2DM51(i),i=1,nlev)
                     Write (unit,101) (CDTM58(i),i=1,nlev)
                     Write (unit,101) (CO2M58(i),i=1,nlev)
                     Write (unit,101) (C2DM58(i),i=1,nlev)

                     Write (unit,102) (STEMP(i),i=1,nlevp1)
                     Write (unit,103) (GTEMP(i),i=1,nlevp1)

                     call close_file (unit)

 101  Format (6F12.8)
 102  Format (5F13.6)
 103  Format (4E16.9)
!-----------------------------------------------------------------------

      End Subroutine Write_CO2_Data

!#######################################################################

      Subroutine Read_CO2_Data (nlev)

      Implicit None
!-----------------------------------------------------------------------
!     Reads co2 transmission functions from file = INPUT/CO2.data
!-----------------------------------------------------------------------
      Integer, Intent(IN) :: nlev
      Integer  unit,nlevp1,i,k
!-----------------------------------------------------------------------

      If (Allocated(CDT51)) Call Error_Mesg ('Read_CO2_Data',  &
                          'CO2 data has already been allocated.', FATAL)

      nlevp1=nlev+1
      unit = open_file (file='INPUT/CO2.data', action='read')

!   B0,B1,B2,B3 ARE COEFFICIENTS USED TO CORRECT FOR THE USE OF 250K IN
!   THE PLANCK FUNCTION USED IN EVALUATING PLANCK-WEIGHTED CO2
!   TRANSMISSION FUNCTIONS. 


      B0 = -.51926410E-4
      B1 = -.18113332E-3
      B2 = -.10680132E-5
      B3 = -.67303519E-7

 101  Format (6f12.8)
 102  Format (5f13.6)
 103  Format (4f16.9)

!  ----- For the 560-800 cm-1 bandwidth -----

      Allocate (CDT51(nlevp1,nlevp1))
      Allocate (CO251(nlevp1,nlevp1))
      Allocate (C2D51(nlevp1,nlevp1))
      Allocate (CDT58(nlevp1,nlevp1))
      Allocate (CO258(nlevp1,nlevp1))
      Allocate (C2D58(nlevp1,nlevp1))

      Do k=1,nlevp1
        Read (unit,101) (CDT51(i,k),i=1,nlevp1)
      EndDo

      Do k=1,nlevp1
        Read (unit,101) (CO251(i,k),i=1,nlevp1)
      EndDo

      Do k=1,nlevp1
        Read (unit,101) (C2D51(i,k),i=1,nlevp1)
      EndDo

      Do k=1,nlevp1
        Read (unit,101) (CDT58(i,k),i=1,nlevp1)
      EndDo

      Do k=1,nlevp1
        Read (unit,101) (CO258(i,k),i=1,nlevp1)
      EndDo

      Do k=1,nlevp1
        Read (unit,101) (C2D58(i,k),i=1,nlevp1)
      EndDo


!  ----- For the 560-670 cm-1 bandwidth -----

      Allocate (CDT31(nlevp1))
      Allocate (CO231(nlevp1))
      Allocate (C2D31(nlevp1))
      Allocate (CDT38(nlevp1))
      Allocate (CO238(nlevp1))
      Allocate (C2D38(nlevp1))

      Read (unit,101) (CDT31(i),i=1,nlevp1)
      Read (unit,101) (CO231(i),i=1,nlevp1)
      Read (unit,101) (C2D31(i),i=1,nlevp1)
      Read (unit,101) (CDT38(i),i=1,nlevp1)
      Read (unit,101) (CO238(i),i=1,nlevp1)
      Read (unit,101) (C2D38(i),i=1,nlevp1)


!  ----- For the 670-800 cm-1 bandwidth -----

      Allocate (CDT71(nlevp1))
      Allocate (CO271(nlevp1))
      Allocate (C2D71(nlevp1))
      Allocate (CDT78(nlevp1))
      Allocate (CO278(nlevp1))
      Allocate (C2D78(nlevp1))

      Read (unit,101) (CDT71(i),i=1,nlevp1)
      Read (unit,101) (CO271(i),i=1,nlevp1)
      Read (unit,101) (C2D71(i),i=1,nlevp1)
      Read (unit,101) (CDT78(i),i=1,nlevp1)
      Read (unit,101) (CO278(i),i=1,nlevp1)
      Read (unit,101) (C2D78(i),i=1,nlevp1)


!  ----- For the 2270-2380 cm-1 bandwidth -----

      Allocate (CO211(nlevp1))
      Allocate (CO218(nlevp1))

      Read (unit,101) (CO211(i),i=1,nlevp1)
      Read (unit,101) (CO218(i),i=1,nlevp1)


! ----- For the 560-800 cm-1 bandwidth -----

      Allocate (CDTM51(nlev))
      Allocate (CO2M51(nlev))
      Allocate (C2DM51(nlev))
      Allocate (CDTM58(nlev))
      Allocate (CO2M58(nlev))
      Allocate (C2DM58(nlev))

      Read (unit,101) (CDTM51(i),i=1,nlev)
      Read (unit,101) (CO2M51(i),i=1,nlev)
      Read (unit,101) (C2DM51(i),i=1,nlev)
      Read (unit,101) (CDTM58(i),i=1,nlev)
      Read (unit,101) (CO2M58(i),i=1,nlev)
      Read (unit,101) (C2DM58(i),i=1,nlev)


!  STEMP IS THE US STANDARD ATMOSPHERES,1976,AT N18 PRESSURES
!  WHERE PSTAR=1013.25 MB

      Allocate (STEMP(nlevp1))

      Read (unit,102) (STEMP(i),i=1,nlevp1)


!  THE WEIGHTING FUNCTION GTEMP=P(K)**0.2*(1.+P(K)/30000.)**0.8/
!  1013250.,WHERE P(K)=PRESSURE,N18 DATA LEVELS FOR PSTAR=
!  1013250.

      Allocate (GTEMP(nlevp1))

      Read (unit,103) (GTEMP(i),i=1,nlevp1)


      call close_file (unit)

!-----------------------------------------------------------------------

      End Subroutine Read_CO2_Data

!#######################################################################

      Subroutine co2_data_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------
      End Subroutine co2_data_init

!#######################################################################

      Subroutine co2_data_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

      End Subroutine co2_data_end

!#######################################################################


      End Module CO2_Data_Mod

