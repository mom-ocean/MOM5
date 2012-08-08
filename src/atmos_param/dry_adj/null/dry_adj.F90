  MODULE DRY_ADJ_MOD

!=======================================================================
!          DRY ADIABATIC ADJUSTMENT       
!=======================================================================

 use       fms_Mod, ONLY: FILE_EXIST, ERROR_MESG, OPEN_NAMELIST_FILE, &
                          CHECK_NML_ERROR, write_version_number, stdlog, &
                          mpp_pe, mpp_root_pe, FATAL, WARNING, CLOSE_FILE
 use Constants_Mod, ONLY: Grav, Kappa
!---------------------------------------------------------------------
 implicit none
 private

 public :: dry_adj, dry_adj_init, dry_adj_end, dry_adj_bdgt

!---------------------------------------------------------------------

 character(len=128) :: version = '$Id: dry_adj.F90,v 10.0 2003/10/24 22:00:29 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized = .false.

!---------------------------------------------------------------------

  contains

!#######################################################################
!#######################################################################

  SUBROUTINE DRY_ADJ ( temp0, pres, pres_int, dtemp, mask )

!=======================================================================
!  DRY ADIABATIC ADJUSTMENT
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     temp0    - Temperature
!     pres     - Pressure
!     pres_int - Pressure at layer interface
!     mask     -  OPTIONAL; floating point mask (0. or 1.) designating 
!                 where data is present
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:,:) :: temp0, pres, pres_int

  real, intent(in), OPTIONAL, dimension(:,:,:) :: mask

!---------------------------------------------------------------------
! Arguments (Intent out)
!     dtemp - Change in temperature
!---------------------------------------------------------------------
  real, intent(out), dimension(:,:,:) :: dtemp

!---------------------------------------------------------------------

      call error_mesg('DRY_ADJ', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE DRY_ADJ

!#####################################################################
!#####################################################################

  SUBROUTINE DRY_ADJ_INIT()

!---------------------------------------------------------------------
! --- WRITE NAMELIST
!---------------------------------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif


!-------------------------------------------------------------------

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('DRY_ADJ_INIT', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE DRY_ADJ_INIT


!#######################################################################
!#######################################################################
  SUBROUTINE DRY_ADJ_END

!-------------------------------------------------------------------

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('DRY_ADJ_END', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE DRY_ADJ_END


!#######################################################################
!#######################################################################

  SUBROUTINE DRY_ADJ_BDGT ( dtemp, pres_int )

!=======================================================================
! Budget check for dry adiabatic adjustment - a debugging tool
!=======================================================================

!---------------------------------------------------------------------
! Arguments (Intent in)
!     dtemp    - Temperature change 
!     pres_int - Pressure at layer interface
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:,:) :: dtemp, pres_int

!---------------------------------------------------------------------

      call error_mesg('DRY_ADJ_BDGT', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE DRY_ADJ_BDGT

!#######################################################################
!#######################################################################
  end MODULE DRY_ADJ_MOD
