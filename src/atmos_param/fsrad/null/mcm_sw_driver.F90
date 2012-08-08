      MODULE MCM_SW_DRIVER_MOD
!
      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

      use mcm_swnew_mod, only: mcm_swnew

implicit none
      private

      character(len=128) :: version = '$Id: mcm_sw_driver.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: siena_201207 $'
      logical            :: module_is_initialized = .false.

      public :: mcm_shortwave_driver, mcm_sw_driver_init, &
                mcm_sw_driver_end

contains

      subroutine mcm_shortwave_driver(                                 &
                     Nclds, KtopSW, KbtmSW, Press, Rh2o, Qo3, CldAmt, &
                     CUVRF, CIRRF, CIRAB, Rco2, CosZ, SSolar,         &
                     Albedo, FSW, DFSW, UFSW, TdtSW, Phalf)


      integer, intent (in), dimension(:,:)     :: Nclds
      integer, intent (in), dimension(:,:,:)      :: KtopSW, KbtmSW
      real, intent (in)   , dimension(:,:,:)      :: Press, Phalf
      real, intent (in)   , dimension(:,:,:)      :: CldAmt, CUVRF,&
                                                  &  CIRRF, CIRAB
      real, intent (in)   , dimension(:,:,:)      :: Rh2o, Qo3

      real, intent (in)                           :: Rco2
      real, intent (in)   , dimension(:,:)        :: CosZ
      real, intent (in)   , dimension(:,:)        :: SSolar
      real, intent (in)   , dimension(:,:)        :: Albedo

      REAL,   INTENT(OUT), DIMENSION(:,:,:)       :: FSW, DFSW, UFSW
      REAL,   INTENT(OUT), DIMENSION(:,:,:)       :: TdtSW


!---------------------------------------------------------------------

      call error_mesg('mcm_shortwave_driver', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_shortwave_driver
! ---------------------------------------------------------------------------------------
      subroutine mcm_sw_driver_init(kx_in)
      integer, intent(in) :: kx_in

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('mcm_sw_driver_init', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_sw_driver_init
! ---------------------------------------------------------------------------------------
      subroutine mcm_sw_driver_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('mcm_sw_driver_end', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_sw_driver_end
! ---------------------------------------------------------------------------------------
      end module MCM_SW_DRIVER_MOD
