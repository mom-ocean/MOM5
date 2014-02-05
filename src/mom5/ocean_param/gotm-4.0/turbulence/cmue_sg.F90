#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The Schumann and Gerz (1995) stability function\label{sec:sg}
!
! !INTERFACE:
   subroutine cmue_sg(nlev)
!
! !DESCRIPTION:
!  This subroutine computes stability functions according to
! \begin{equation}
! c_{\mu}=c_{\mu}^0,\qquad c'_{\mu}=\frac{c_{\mu}^0}{Pr_t}
! \end{equation}
! with constant $c_{\mu}^0$. Based simulation data on stratified homogeneous
! shear-flows, \cite{SchumannGerz95} proposed the empirical relation
! for the turbulent Prandtl--number,
! \begin{equation}
!   Pr_t = Pr_t^0 \exp\left(-\frac{Ri}{Pr_t^0 Ri^{\infty}}\right)
!   -\frac{Ri}{Ri^{\infty}}
!   \comma
! \end{equation}
! where where $Ri$ is the gradient Richardson--number and $Pr_t^0$
! is the turbulent Prandtl--number for $Ri \rightarrow 0$. $Pr_t^0$
! and the fixed value $c_\mu^0$ have to be set in {\tt gotmturb.nml}.
! \cite{SchumannGerz95}  suggested $Pr_t^0=0.74$ and $Ri^{\infty}=0.25$.
!
! !USES:
   use turbulence, only: Prandtl0_fix,cm0_fix
   use turbulence, only: cmue1,cmue2,as,an
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: cmue_sg.F90,v $
!  Revision 20.0  2013/12/14 00:13:34  fms
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
!  Revision 1.9  2010-09-17 12:53:52  jorn
!  extensive code clean-up to ensure proper initialization and clean-up of all variables
!
!  Revision 1.8  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.7  2005/11/15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.6  2005/07/18 08:54:56  lars
!  changed docu for html compliance
!
!  Revision 1.5  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2004/08/18 12:53:07  lars
!  updated documentation
!
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 09:02:04  gotm
!  Added new Generic Turbulence Model +
!  improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: Ri,Prandtl
   REALTYPE,parameter        :: limit=3.
!
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      Ri=an(i)/(as(i)+1.e-8)   ! Gradient Richardson number
      if (Ri.ge.1e-10) then
         Prandtl=Prandtl0_fix*exp(-Ri/(Prandtl0_fix*0.25))+Ri/0.25
      else
         Prandtl=Prandtl0_fix
      end if

      cmue1(i)=cm0_fix
      cmue2(i)=cm0_fix/min(limit,Prandtl)

   end do
   return
   end subroutine cmue_sg
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
