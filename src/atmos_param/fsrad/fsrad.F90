                         Module FSrad_Mod

!-----------------------------------------------------------------------
!-------------------- PUBLIC Radiation routines ------------------------

      Use        MCM_LW_Mod, ONLY: MCM_LW_Rad
      Use MCM_SW_Driver_Mod, ONLY: mcm_shortwave_driver

      Use   ShortWave_Mod, ONLY: SWRad
      Use    LongWave_Mod, ONLY: LWRad, Rad_DeAlloc
      Use      RdParm_Mod, ONLY: RdParm_Init
      Use    Rad_Diag_Mod, ONLY: Radiag
      Use    CO2_Data_Mod, ONLY: CO2_Data

      Use         Fms_Mod, ONLY: mpp_pe, mpp_root_pe, write_version_number, &
                                 error_mesg, FATAL
      Use   Constants_Mod, ONLY: stefan

      implicit none
      private

      public  FSrad, RdParm_Init, CO2_Data, fsrad_init, fsrad_end
!-----------------------------------------------------------------------

      character(len=128) :: version = '$Id: fsrad.F90,v 14.0 2007/03/15 22:03:19 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical            :: module_is_initialized = .false.

      real, parameter :: Day_Length=86400.
      real, parameter :: RATco2MW=1.519449738
      real, parameter :: PerSec=1./Day_Length

!-----------------------------------------------------------------------

CONTAINS

!#######################################################################

      Subroutine FSrad (ip,jp,Press,Temp,Rh2o,Qo3,              &
                        phalf,do_mcm_radiation,             &
                        Nclds,KtopSW,KbtmSW,Ktop,Kbtm,CldAmt,   &
                        EmCld,CUVRF,CIRRF,CIRAB,Albedo,RVco2,   &
                        CosZ,Solar,                             &
                        SWin,SWout,OLR,SWupS,SWdnS,LWupS,LWdnS, &
                        TdtSW,TdtLW, Ksfc,Psfc)

!-----------------------------------------------------------------------
Integer, Intent(IN)                    :: ip,jp
   Real, Intent(IN), Dimension(:,:,:)  :: Press,Temp,Rh2o,Qo3
   Real, Intent(IN), Dimension(:,:,:)  :: phalf
Logical, Intent(IN)                    :: do_mcm_radiation
Integer, Intent(IN), Dimension(:,:)    :: Nclds
Integer, Intent(IN), Dimension(:,:,:)  :: KtopSW,KbtmSW,Ktop,Kbtm
   Real, Intent(IN), Dimension(:,:,:)  :: CldAmt,EmCld,CUVRF,CIRRF,CIRAB
   Real, Intent(IN), Dimension(:,:)    :: Albedo,CosZ,Solar
   Real, Intent(IN)                    :: RVco2

   Real, Intent(OUT), Dimension(:,:)   :: SWin,SWout,OLR,SWupS,SWdnS,  &
                                                       LWupS,LWdnS
   Real, Intent(OUT), Dimension(:,:,:) :: TdtSW,TdtLW

Integer, Intent(IN), Dimension(:,:), Optional :: Ksfc
   Real, Intent(IN), Dimension(:,:), Optional :: Psfc
!-----------------------------------------------------------------------
   Real, Dimension(Size(Rh2o,1),Size(Rh2o,2),Size(Rh2o,3)+1) ::  &
                                    FSW,DFSW,UFSW
   Real, Dimension(Size(Rh2o,1),Size(Rh2o,2)) :: SSolar,GrnFlux,TopFlux
   Real  Rco2
Logical  SunUp
Integer  i,j,IX,JX,KX
!-----------------------------------------------------------------------


      IX=Size(Rh2o,1)
      JX=Size(Rh2o,2)
      KX=Size(Rh2o,3)

      SunUp=.false.
      Do j=1,JX
      Do i=1,IX
         If (CosZ(i,j) > 0.0) Then
            SunUp=.true.
            EXIT
         EndIf
      EndDo
      EndDo

!-----------------------------------------------------------------------
!----- convert solar constant from W/m2 to ly/min ------

      SSolar = Solar * (6./4186.)

!-----------------------------------------------------------------------
!----------------------- Shortwave Radiation ---------------------------

      If (SunUp) Then
         Rco2=RVco2*RATco2MW

        if ( do_mcm_radiation ) then
         Call mcm_shortwave_driver &
                    (Nclds, KtopSW, KbtmSW, Press, Rh2o, Qo3, CldAmt, &
                     CUVRF, CIRRF, CIRAB, Rco2, CosZ, SSolar,         &
                     Albedo, FSW, DFSW, UFSW, TdtSW, phalf )
        else
         Call SWRad (Nclds, KtopSW, KbtmSW, Press, Rh2o, Qo3, CldAmt, &
                     CUVRF, CIRRF, CIRAB, Rco2, CosZ, SSolar,         &
                     Albedo, FSW, DFSW, UFSW, TdtSW, Ksfc, Psfc       )
        endif

      Else
         FSW=0.0; DFSW=0.0; UFSW=0.0; TdtSW=0.0
      EndIf


!-----------------------------------------------------------------------
!----------------------- Longwave Radiation ----------------------------

      if ( do_mcm_radiation ) then
        Call MCM_LW_Rad (Ktop, Kbtm, Nclds, EmCld, Press, Temp, Rh2o, Qo3,&
                    CldAmt, RVco2, TdtLW, GrnFlux, TopFlux, phalf )
      else
        Call LWRad (Ktop, Kbtm, Nclds, EmCld, Press, Temp, Rh2o, Qo3,  &
                    CldAmt, RVco2, TdtLW, GrnFlux, TopFlux, Ksfc, Psfc )
      endif

!-----------------------------------------------------------------------
!----------------------- Radiation Diagnostics -------------------------
      If (ip > 0 .and. jp > 0) Then
         Call Radiag (Press,Temp,Rh2o,Rco2,Qo3,CldAmt,Ktop,Kbtm,Nclds, &
                      TdtLW,GrnFlux,  FSW,DFSW,UFSW,TdtSW,  &
                      KtopSW,KbtmSW,EmCld,CUVRF,CIRRF,CIRAB,  &
                      Albedo,CosZ,SSolar,   ip,jp)
      EndIf
!-----------------------------------------------------------------------
      if(.not.do_mcm_radiation ) then
        Call Rad_DeAlloc
      endif
!-----------------------------------------------------------------------
!    **** Output fluxes at the top and bottom of atmosphere ****
!            **** convert ergs/cm2/s to watts/m2 ****

!  --- TOA ---
          SWin (:,:)=DFSW(:,:,1)  *1.E-3
          SWout(:,:)=UFSW(:,:,1)  *1.E-3
          OLR  (:,:)=TopFlux(:,:) *1.E-3

!  --- SFC ---
   If (Present(Ksfc) .and. Present(Psfc)) Then
       Do j=1,JX
       Do i=1,IX
          SWupS(i,j)=UFSW(i,j,Ksfc(i,j)+1)*1.E-3
          SWdnS(i,j)=DFSW(i,j,Ksfc(i,j)+1)*1.E-3
          LWupS(i,j)=stefan*Temp(i,j,Ksfc(i,j)+1)**4
       EndDo
       EndDo
   Else
          SWupS(:,:)=UFSW(:,:,KX+1)*1.E-3
          SWdnS(:,:)=DFSW(:,:,KX+1)*1.E-3
          LWupS(:,:)=stefan*Temp(:,:,KX+1)**4
   EndIf
          LWdnS(:,:)=LWupS(:,:)-GrnFlux(:,:)*1.E-3


!-----------------------------------------------------------------------

      TdtSW=TdtSW*PerSec
      TdtLW=TdtLW*PerSec

!-----------------------------------------------------------------------

      End Subroutine FSrad

!#######################################################################

      Subroutine fsrad_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      End Subroutine fsrad_init

!#######################################################################

      Subroutine fsrad_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      End Subroutine fsrad_end

!#######################################################################



                     End Module FSrad_Mod

