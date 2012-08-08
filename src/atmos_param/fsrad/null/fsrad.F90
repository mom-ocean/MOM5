
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

      character(len=128) :: version = '$Id: fsrad.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: siena_201207 $'
      logical            :: module_is_initialized = .false.

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
Integer, Intent(IN)                     :: ip,jp
   Real, Intent(IN),  Dimension(:,:,:)  :: Press,Temp,Rh2o,Qo3
   Real, Intent(IN),  Dimension(:,:,:)  :: phalf
Logical, Intent(IN)                     :: do_mcm_radiation
Integer, Intent(IN),  Dimension(:,:)    :: Nclds
Integer, Intent(IN),  Dimension(:,:,:)  :: KtopSW,KbtmSW,Ktop,Kbtm
   Real, Intent(IN),  Dimension(:,:,:)  :: CldAmt,EmCld,CUVRF,CIRRF,CIRAB
   Real, Intent(IN),  Dimension(:,:)    :: Albedo,CosZ,Solar
   Real, Intent(IN)                     :: RVco2

   Real, Intent(OUT), Dimension(:,:)    :: SWin,SWout,OLR,SWupS,SWdnS,  &
                                                       LWupS,LWdnS
   Real, Intent(OUT), Dimension(:,:,:)  :: TdtSW,TdtLW

Integer, Intent(IN),  Dimension(:,:), Optional :: Ksfc
   Real, Intent(IN),  Dimension(:,:), Optional :: Psfc
!---------------------------------------------------------------------

      call error_mesg('FSrad', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine FSrad

!#######################################################################

      Subroutine fsrad_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('fsrad_init', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine fsrad_init

!#######################################################################

      Subroutine fsrad_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('fsrad_end', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine fsrad_end

!#######################################################################


                     End Module FSrad_Mod

