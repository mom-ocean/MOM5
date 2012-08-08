      MODULE MCM_LW_MOD

!   Added interface routine (lw_rad_ss) which is called by
!     fsrad and which calls lwcool in this
!     module after constructing the appropriate inputs.

      USE Constants_Mod, ONLY: grav, tfreeze

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

implicit none
      private

!------------ VERSION NUMBER ----------------

      character(len=128) :: version = '$Id: mcm_lw.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: siena_201207 $'
      logical            :: module_is_initialized = .false.

      public  MCM_LW_RAD, mcm_lw_init, mcm_lw_end

!     -------------------------------------------------

!-----------------------------------------------------------------------
!--------------------- G L O B A L   D A T A ---------------------------
!-----------------------------------------------------------------------


      contains

!#######################################################################
      subroutine mcm_lw_init(ix_in, jx_in, kx_in)
      integer, intent(in) :: ix_in, jx_in, kx_in

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('mcm_lw_init', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_lw_init
!#######################################################################

      subroutine MCM_LW_END

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('MCM_LW_END', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine MCM_LW_END

!#######################################################################
      SUBROUTINE MCM_LW_RAD (KTOP,KBTM,NCLDS,EMCLD, &
                      PRES,TEMP,RH2O,QO3,CAMT, &
                      RRVCO2,  HEATRA,GRNFLX,TOPFLX, phalf)

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)    :: NCLDS
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRES,TEMP
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: CAMT
      REAL,    INTENT(IN)                      :: RRVCO2
 
      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: HEATRA
      REAL,   INTENT(OUT), DIMENSION(:,:)    :: GRNFLX,TOPFLX

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: phalf



!---------------------------------------------------------------------

      call error_mesg('MCM_LW_RAD', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE MCM_LW_RAD

      end module mcm_lw_mod
