      module mcm_swnew_mod

      use mcm_swtbls_mod, only: aaa, aab
      Use       Fms_Mod, ONLY: write_version_number, mpp_pe, mpp_root_pe, &
                               error_mesg, FATAL

implicit none 
private 

      character(len=128) :: version = '$Id: mcm_swnew.F90,v 11.0 2004/09/28 19:18:44 fms Exp $'
      character(len=128) :: tagname = '$Name: siena_201207 $'
      logical            :: module_is_initialized = .false.

public mcm_swnew, mcm_swnew_init, mcm_swnew_end

contains

      subroutine mcm_swnew( cosz, rco2, rh2o, ro3, pp, &
              cwca, cwcb, coca, cloudy, kthsw, kbhsw, ssolar, pr2, &
              flx, heat, grdflx, ncv, kx, UF, DF)

      integer, intent (in)                       :: ncv, kx
      real   , intent (in)  :: cosz, rco2, Ssolar

      real   , intent (in), dimension(kx)        :: rh2o, ro3
      real   , intent (in), dimension(0:kx)      :: pp
      real   , intent (in), dimension(1:kx+2)    :: cwca, cwcb, coca, cloudy
      real   , intent (in), dimension(1:kx+1)    :: pr2

      integer, intent (in), dimension(1:kx+2)    :: kthsw, kbhsw

      real   , intent (out)                      :: grdflx
      real   , intent (out), dimension(1:kx+1)   :: flx, heat, UF, DF

!---------------------------------------------------------------------

      call error_mesg('mcm_swnew', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_swnew
! ---------------------------------------------------------------------------------------
      subroutine mcm_swnew_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('mcm_swnew_init', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_swnew_init
! ---------------------------------------------------------------------------------------
      subroutine mcm_swnew_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('mcm_swnew_end', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_swnew_end

      end module mcm_swnew_mod
