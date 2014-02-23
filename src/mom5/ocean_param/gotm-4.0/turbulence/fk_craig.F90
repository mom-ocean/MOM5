#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: TKE flux from wave-breaking\label{sec:fkCraig}
!
! !INTERFACE:
   REALTYPE  function fk_craig(u_tau,eta)
!
! !DESCRIPTION:
! This functions returns the flux of $k$ caused by breaking surface waves
! according to
! \begin{equation}
!  \label{craig}
!   F_k = \eta u_*^3
!  \point
! \end{equation}
! This form has also been used by \cite{CraigBanner94}, who suggested
! $\eta \approx 100$.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: u_tau
   REALTYPE, intent(in)                :: eta
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: fk_craig.F90,v $
!  Revision 20.0  2013/12/14 00:13:39  fms
!  Merged revision 1.1.2.1 onto trunk
!
!  Revision 1.1.2.1  2012/05/15 16:00:53  smg
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
!  Revision 1.5  2007-07-23 11:28:39  hb
!  cw for Craig-Banner wave breaking from namelist now used in fk_craig.F90
!
!  Revision 1.4  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2004/08/18 12:50:57  lars
!  updated documentation
!
!  Revision 1.2  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.1  2003/03/10 09:00:36  gotm
!  Part of new generic turbulence model
!
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fk_craig = eta*u_tau**3.

   end function fk_craig
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
