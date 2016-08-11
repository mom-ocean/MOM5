MODULE mg_const_mod

use  fms_mod,            only :  write_version_number
use  constants_mod,      only :  pi

implicit none
private

!-------------------------------------------------------------------------
!----interfaces-----------------------------------------------------------

public mg_const_init

!------------------------------------------------------------------------
!       DECLARE VERSION NUMBER
!------------------------------------------------------------------------
Character(len=128) :: Version = '$Id: mg_const.F90,v 20.0 2013/12/13 23:21:55 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'
 
!-------------------------------------------------------------------------
!---module variables------------------------------------------------------

INTEGER, PUBLIC, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)

!double precision
INTEGER, PUBLIC, PARAMETER :: mg_pr=dp


! bulk density ice (kg/m3) from Reisner et al. (1998)
REAL(kind=mg_pr), PUBLIC, PARAMETER :: rhoi = 500._mg_pr    
! bulk density liquid (kg/m3) from Reisner et al. (1998) 
REAL(kind=mg_pr), PUBLIC, PARAMETER :: rhow = 1000._mg_pr   

! cloud ice mass-diameter relationship

REAL(kind=mg_pr), PUBLIC, PARAMETER :: ci = rhoi*pi/6._mg_pr
REAL(kind=mg_pr), PUBLIC, PARAMETER :: di = 3._mg_pr

REAL(kind=mg_pr), PUBLIC, PARAMETER :: di_mg = di
REAL(kind=mg_pr), PUBLIC, PARAMETER :: ci_mg = ci





logical   :: module_is_initialized = .false.


CONTAINS

!#########################################################################

subroutine mg_const_init

      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    write version number to output file.
!-------------------------------------------------------------------------
      call write_version_number(version, tagname)

!-------------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------------


end subroutine mg_const_init

!#########################################################################


END MODULE mg_const_mod
