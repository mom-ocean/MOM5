
                        MODULE SHORTWAVE_MOD

!-----------------------------------------------------------------------

      USE  RDPARM_MOD, ONLY: LMAX,LP1,LLP1,LP2,LLP2,NB
      USE  HCONST_MOD, ONLY: DIFFCTR,GINV,O3DIFCTR,RADCON

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

implicit none
private

!------- interfaces -------
      PUBLIC  SWRAD, SHORTWAVE_INIT, SHORTWAVE_END

      character(len=128) :: version = '$Id: shortwave.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: siena_201207 $'
      logical            :: module_is_initialized = .false.

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

!---------------------------------------------------------------------

      call error_mesg('SWRAD', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE SWRAD
      
!#######################################################################

      SUBROUTINE SHORTWAVE_INIT

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('SHORTWAVE_INIT', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE SHORTWAVE_INIT

!#######################################################################

      SUBROUTINE SHORTWAVE_END

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('SHORTWAVE_END', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE SHORTWAVE_END

!#######################################################################

      END MODULE SHORTWAVE_MOD
