#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update dimensionless alpha's\label{sec:alpha}
!
! !INTERFACE:
   subroutine alpha_mnb(nlev,NN,SS)
!
! !DESCRIPTION:
! This subroutine updates the dimensionless numbers $\alpha_M$, $\alpha_N$,
! and $\alpha_b$ according to \eq{alphaMN}. Note that according to \eq{Nbar}
! and \eq{NbarVertical} the following identities are valid
! \begin{equation}
!  \label{alphaIdentities}
!    \alpha_M = \overline{S}^2 \comma
!    \alpha_N = \overline{N}^2 \comma
!    \alpha_b = \overline{T}   \point
! \end{equation}
!
!
! !USES:
  use turbulence,  only:     tke,eps,kb
  use turbulence,  only:     as,an,at
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)      :: nlev
  REALTYPE, intent(in)      :: NN(0:nlev),SS(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: alpha_mnb.F90,v $
!  Revision 20.0  2013/12/14 00:13:27  fms
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
!  Revision 1.2  2006-03-20 09:06:37  kbk
!  removed explicit double precission dependency
!
!  Revision 1.1  2005/06/27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer              :: i
  REALTYPE             :: tau2

!-----------------------------------------------------------------------
!BOC

  do i=0,nlev
     tau2   = tke(i)*tke(i) / ( eps(i)*eps(i) )
     as(i)  = tau2 * SS(i)
     an(i)  = tau2 * NN(i)
     at(i)  = tke(i)/eps(i) * kb(i)/eps(i)

!    clip negative values
     as(i) = max(as(i),1.e-10*_ONE_)
     at(i) = max(at(i),1.e-10*_ONE_)
  end do

  return
end subroutine alpha_mnb

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
