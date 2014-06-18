#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic kb-equation\label{sec:kbalgebraic}
!
! !INTERFACE:
   subroutine kbalgebraic(nlev)
!
! !DESCRIPTION:
! The algebraic equation for $k_b$ simply assumes equilibrium in \eq{kbeq},
! \begin{equation}
!   \label{kbEquilibrium}
!   P_b = \epsilon_b
!   \point
! \end{equation}
! This equation can be re-written as
! \begin{equation}
!   \label{kbAgebraic}
!   k_b = \dfrac{k_b \epsilon}{k \epsilon_b} \dfrac{k}{\epsilon} P_b
!       = r \dfrac{k}{\epsilon} P_b = c_b \dfrac{k}{\epsilon} P_b
!   \comma
! \end{equation}
! where we used the definition of the time scale ratio $r$ in
! \eq{DefR}, and assumed that $r=c_b$ is a constant.
!

!
! !USES:
   use turbulence,  only:     tke,eps,kb,Pb
   use turbulence,  only:     ctt,kb_min

  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
   integer,  intent(in)                 :: nlev

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: kbalgebraic.F90,v $
!  Revision 20.0  2013/12/14 00:13:44  fms
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
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

   integer                             :: i

!-----------------------------------------------------------------------
!BOC

   do i=0,nlev
      kb(i) = ctt*tke(i)/eps(i)*Pb(i)

      !  clip at kb_min
      kb(i) = max(kb(i),kb_min)
   enddo

   return
   end subroutine kbalgebraic

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
