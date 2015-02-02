#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update time scale ratio
!
! !INTERFACE:
   subroutine r_ratio(nlev)
!
! !DESCRIPTION:
! This routine updates the ratio $r$ of the dissipation
! time scales as defined in \eq{DefR}.
!
! !USES:
  use turbulence,  only:     tke,eps,kb,epsb
  use turbulence,  only:     r

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in)        :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: r_ratio.F90,v $
!  Revision 20.0  2013/12/14 00:13:51  fms
!  Merged revision 1.1.2.1 onto trunk
!
!  Revision 1.1.2.1  2012/05/15 16:00:54  smg
!  initial cvs ci for these modules to mom5.
!  AUTHOR:Griffies
!  REVIEWERS:
!  TEST STATUS:
!  CHANGES PUBLIC INTERFACES?
!  CHANGES ANSWERS?
!
!  Revision 1.1.2.1.390.1  2012/04/23 20:30:29  smg
!  updated to the gotm-2012.03.09 CVS tag.
!  AUTHOR:Martin Schmidt
!  REVIEWERS:Griffies
!  TEST STATUS:
!  CHANGES PUBLIC INTERFACES?
!  CHANGES ANSWERS?
!
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
!BOC

!KBK   r = kb/epsb*eps/tke
   r = kb*eps/(epsb*tke)

   return
end subroutine r_ratio

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
